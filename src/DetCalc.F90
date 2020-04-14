#include "macros.h"
MODULE DetCalc
        use constants, only: dp,n_int
        use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB, tStoreSpinOrbs, &
             tGAS, t_non_hermitian
        use disconnected_gasci, only: init_disconnected_GAS
        use gasci, only: init_GAS, GAS_specification, GAS_exc_gen, possible_GAS_exc_gen, operator(==)
        use sort_mod

        use DetCalcData

        use MemoryManager, only: TagIntType

        use gndts_mod, only: gndts

        use UMatCache, only: UMat2D, tUMat2D, tDeferred_UMat2D, SetupUMat2d_dense

        use procedure_pointers, only: get_umat_el

    IMPLICIT NONE
     save

!From input
      INTEGER DETINV !The index in the list of dets of a det to investigate
      INTEGER IOBS,JOBS,KOBS

      LOGICAL TCALCHMAT,TENERGY,TREAD,TBLOCK
      LOGICAL tFindDets           !Set if we are to enumerate all determinants within given constraints
      LOGICAL tCompressDets       !Set if once we've found the dets we compress to bit format

      TYPE(BasisFN), pointer :: BLOCKSYM(:)  !The Symmetry of each block.  nBlocks elements
      INTEGER(TagIntType) :: tagBlockSym
      INTEGER,ALLOCATABLE :: NBLOCKSTARTS(:) !Index of the first det of different symmetry blocks in the complete list of dets
      INTEGER(TagIntType) :: tagNBLOCKSTARTS=0
      INTEGER NBLOCKS                        !Number of Symmetry blocks
      INTEGER iFDet                       ! The index of the Fermi det in the list of dets.
      HElement_t(dp), pointer :: CKN(:,:) !  (nDet,nEval)  Temporary storage for the Lanczos routine
      INTEGER(TagIntType) :: tagCKN=0

      real(dp) , ALLOCATABLE :: ExpandedHamil(:,:)    ! (NDet,NDet) This is the hamiltonian in expanded form,
                                                      !so that it can be histogrammed against.

      INTEGER iExcitLevel                 ! The excitation level at which determinants are cut off.

CONTAINS
    Subroutine DetCalcInit

        Use global_utilities
        Use Determinants, only:  FDet, specdet, tSpecDet, tDefineDet, &
                                 DefDet, write_det
        Use IntegralsData, only : NFROZEN
        use SystemData, only : tCSFOLD,lms, lms2, nBasis, nBasisMax, nEl, SymRestrict
        use SystemData, only : Alat, arr, brr, boa, box, coa, ecore, g1,Beta
        use SystemData, only : tParity, tSpn,Symmetry,STot, NullBasisFn, tUHF,tMolpro
        use sym_mod
        use LoggingData,    only : tLogDets
        use legacy_data, only: irat
        use HElem
        use util_mod, only: get_free_unit, NECI_ICOPY
        Type(BasisFn) ISym

        integer i, j, ii, iunit
        integer ierr, norb
        integer nDetTot
        logical isvaliddet

        character(25), parameter :: this_routine='DetCalcInit'


      IF(.NOT.TCALCHMAT) THEN
         WRITE(6,*) "Not storing the H matrix."
         IF(TENERGY.AND..NOT.TBLOCK) THEN
            WRITE(6,*) "Cannot calculate energies without blocking the Hamiltonian."
            TENERGY=.FALSE.
         ENDIF
         IF(TENERGY.AND.NBLK.NE.0) THEN
!C.. We're doing a Lanczos without calculating the H mat
            WRITE(6,*) "Cannot perform Lanczos without Hamiltonian"
            TENERGY=.FALSE.
         ENDIF
      ENDIF

      ! If we want to have UMat2D, and it isn't yet filled in, generate it
      ! here. All of the integrals setup/freezing etc is done...
      if (tDeferred_Umat2d .and. .not. tUMat2D) then

          ASSERT(.not. btest(nbasis, 0))

          ! And fill in the array
          call SetupUMat2d_dense(nBasis)
      end if


!Copied Specdet information from Calc.F, so if inspect is present, but no determinant/csf specified, it will still run.
      if(TSPECDET) then
         if(TCSFOLD) then
             WRITE(6,*) "TSPECDET set.  SPECDET is"
             call write_det (6, SPECDET, .true.)
             CALL NECI_ICOPY(NEL,SPECDET,1,FDET,1)
             CALL GETCSFFROMDET(FDET,SPECDET,NEL,STOT,LMS)
             WRITE(6,*) "CSF with 2S=",STOT," and 2Sz=",LMS," now in SPECDET is"
             call write_det (6, SPECDET, .true.)
         elseif(.not.associated(specdet)) then
            !specdet not allocated. Allocate it and copy fdet
             allocate(specdet(nel))
             WRITE(6,*) "TSPECDET set, but not allocated.  using FDET"
             CALL NECI_ICOPY(NEL,FDET,1,SPECDET,1)
         elseif(.not.ISVALIDDET(SPECDET,NEL)) then
             WRITE(6,*) "TSPECDET set, but invalid.  using FDET"
!             tSpecDet=.false.
             CALL NECI_ICOPY(NEL,FDET,1,SPECDET,1)
         endif
      ELSEIF(TCSFOLD) THEN  !No help given on generating this CSF.  Let's just get a single one out of GNCSFs
         NDET=1
         CALL GNCSFS(NEL,nBasis,BRR,NBASISMAX,FDET,.FALSE.,G1,TSPN,LMS2,TPARITY, &
     &      SymRestrict,NDET,IFDET,.FALSE.,0,II,.FALSE.,0)   !II is just a dummmy to receive IC
         WRITE(6,*) "CSF with 2S=",STOT," and 2Sz=",LMS," now in FDET is"
         IFDET=0
         NDET=0
         call write_det (6, FDET, .true.)
      ENDIF


!C      IF(TCALCHMAT.OR.NPATHS.NE.0.OR.DETINV.GT.0.OR.TBLOCK) THEN
      iExcitLevel=ICILEVEL
      IF(tFindDets) THEN
!C..Need to determine the determinants
         IF(iExcitLevel.NE.0) THEN
            WRITE(6,*) "Performing truncated CI at level ",iExcitLevel
            IF(TSPECDET) THEN
               WRITE(6,*) "Using SPECDET:"
               call write_det (6, SPECDET, .true.)!
               CALL NECI_ICOPY(NEL,SPECDET,1,FDET,1)
            ELSE
               WRITE(6,*) "Using Fermi DET:"
               call write_det (6, FDET, .true.)
            ENDIF
            IF(TCSFOLD) WRITE(6,*) "Determining CSFs."
!C.. if we're doing a truncated CI expansion
            CALL GENEXCIT(FDET,iExcitLevel,NBASIS,NEL,0,(/0.0_dp/),NDET,1,G1,.TRUE.,NBASISMAX,.TRUE.)
            WRITE(6,*) "NDET out of GENEXCIT ",NDET
!C.. We need to add in the FDET
            NDET=NDET+1
            II=NDET
            NBLOCKS=1
         ELSEIF(TBLOCK) THEN
            WRITE(6,*) "Determining determinants and blocks."
            IF(TPARITY) THEN
               WRITE(6,*) "Using symmetry restriction:"
               CALL WRITEALLSYM(6,SymRestrict,.TRUE.)
            ENDIF
            IF(TSPN) THEN
               WRITE(6,*) "Using spin restriction:",LMS
            ENDIF
            if(tUHF.and.tMolpro) then
                !When breaking spin symmetry in molpro, it is important to occupy alpha orbs preferentially
                CALL GNDTS_BLK(NEL,nBasis,BRR,NBASISMAX,NMRKS, .TRUE.,             &
     &                NDET,G1,II,NBLOCKSTARTS,NBLOCKS,TSPN,-LMS2,TPARITY,        &
     &               SymRestrict,IFDET,.NOT.TREAD,NDETTOT,BLOCKSYM)
            else
                CALL GNDTS_BLK(NEL,nBasis,BRR,NBASISMAX,NMRKS, .TRUE.,             &
     &                NDET,G1,II,NBLOCKSTARTS,NBLOCKS,TSPN,LMS2,TPARITY,        &
     &               SymRestrict,IFDET,.NOT.TREAD,NDETTOT,BLOCKSYM)
            endif
            WRITE(6,*) "NBLOCKS:",NBLOCKS
         ELSEIF(TCSFOLD) THEN
            WRITE(6,*) "Determining CSFs."
!C determinants.
            WRITE(6,*) "Determining determinants and blocks."
            IF(TPARITY) THEN
               WRITE(6,*) "Using symmetry restriction:"
               CALL WRITEALLSYM(6,SymRestrict,.TRUE.)
            ENDIF
            IF(TSPN) THEN
               WRITE(6,*) "Using spin restriction:",LMS
            ENDIF
            NDET=0
            CALL GNCSFS(NEL,nBasis,BRR,NBASISMAX,NMRKS,.TRUE.,G1,TSPN,LMS2,TPARITY,        &
     &         SymRestrict,NDET,IFDET,.FALSE.,0,0,.FALSE.,0)
            NBLOCKS=1
            II=NDET
         ELSE
            WRITE(6,*) "Determining determinants."
            IF(TPARITY) THEN
               WRITE(6,*) "Using symmetry restriction:"
               CALL WRITEALLSYM(6,SymRestrict,.TRUE.)
            ENDIF
            IF(TSPN) THEN
               WRITE(6,*) "Using spin restriction:",LMS
            ENDIF
            if(tUHF.and.tMolpro) then
                !When breaking spin symmetry in molpro, it is important to occupy alpha orbs preferentially
                CALL GNDTS(NEL,nBasis,BRR,NBASISMAX,NMRKS,.TRUE.,G1,TSPN,-LMS,TPARITY,SymRestrict,II,IFDET)
            else
                CALL GNDTS(NEL,nBasis,BRR,NBASISMAX,NMRKS,.TRUE.,G1,TSPN,LMS,TPARITY,SymRestrict,II,IFDET)
            endif
            NBLOCKS=1
            NDET=II
         ENDIF
!C..
         IF(II.EQ.0) THEN
            WRITE(6,*) "No determinants found.  Cannot continue"
            call stop_all(this_routine, "No determinants found.  Cannot continue")
         ENDIF
!C.. NEL now only includes active electrons
         WRITE(6,*) "Number of determinants found to be: ",II
         WRITE(6,*) "Allocating initial memory for calculation of energy..."
         CALL neci_flush(6)
         Allocate(NMrks(nEl,II),stat=ierr)
         LogAlloc(ierr,'NMRKS',NEL*II,4,tagNMRKS)
         NMRKS(1:NEL,1:II)=0
         allocate(NBLOCKSTARTS(NBLOCKS+1),stat=ierr)
         call LogMemAlloc('NBLOCKSTARTS',NBLOCKS+1,4,this_routine,tagNBLOCKSTARTS,ierr)
         NBLOCKSTARTS(1:NBLOCKS+1)=0
         Allocate(BlockSym(NBLOCKS+1),stat=ierr)
         LogAlloc(ierr, 'BLOCKSYM', NBLOCKS+1,BasisFNSizeB, tagBlockSym)

         BLOCKSYM(1:NBLOCKS)=NullBasisFn
!C..



         NDET=II
         IF(iExcitLevel.NE.0) THEN
!C.. Use HAMIL to temporarily hold a list of excitation levels
            CALL NECI_ICOPY(NEL,FDET,1,NMRKS,1)
            Allocate(Hamil(II), stat=ierr)
            LogAlloc(ierr, 'HAMIL', II, HElement_t_sizeB, tagHamil)
            NDET=0
            CALL GENEXCIT(FDET,iExcitLevel,NBASIS,NEL,NMRKS(1,2),HAMIL,NDET,1,G1,.TRUE.,NBASISMAX,.FALSE.)
            Deallocate(Hamil)
            LogDealloc(tagHamil)
            NDET=NDET+1
            NBLOCKSTARTS(1)=1
            NBLOCKSTARTS(2)=II+1
            IFDET=1
         ELSEIF(TBLOCK) THEN
            if(tUHF.and.tMolpro) then
                !When breaking spin symmetry in molpro, it is important to occupy alpha orbs preferentially
                CALL GNDTS_BLK(NEL,nBasis,BRR,NBASISMAX,NMRKS, .FALSE.,NDET,G1,II,NBLOCKSTARTS,NBLOCKS,TSPN,-LMS2,TPARITY, &
     &               SymRestrict,IFDET,.NOT.TREAD,NDETTOT,BLOCKSYM)
            else
                CALL GNDTS_BLK(NEL,nBasis,BRR,NBASISMAX,NMRKS, .FALSE.,NDET,G1,II,NBLOCKSTARTS,NBLOCKS,TSPN,LMS2,TPARITY, &
     &               SymRestrict,IFDET,.NOT.TREAD,NDETTOT,BLOCKSYM)
            endif
         ELSEIF(TCSFOLD) THEN
            NDET=0  !This will be reset by GNCSFS
            CALL GNCSFS(NEL,nBasis,BRR,NBASISMAX,NMRKS,.FALSE.,G1,TSPN,LMS2,TPARITY, &
     &         SymRestrict,NDET,IFDET,.FALSE.,0,0,.FALSE.,0)
               NBLOCKSTARTS(1)=1
               NBLOCKSTARTS(2)=II+1
         ELSE
            if(tUHF.and.tMolpro) then
                !When breaking spin symmetry in molpro, it is important to occupy alpha orbs preferentially
                CALL GNDTS(NEL,nBasis,BRR,NBASISMAX,NMRKS,.FALSE.,G1,TSPN,-LMS,TPARITY,SymRestrict,II,IFDET)
            else
                CALL GNDTS(NEL,nBasis,BRR,NBASISMAX,NMRKS,.FALSE.,G1,TSPN,LMS,TPARITY,SymRestrict,II,IFDET)
            endif
               NBLOCKSTARTS(1)=1
               NBLOCKSTARTS(2)=II+1
         ENDIF
         if(tLogDets) THEN
            iunit = get_free_unit()
            OPEN(iunit,FILE='DETS',STATUS='UNKNOWN')
            DO I=1,NDET
               call write_det (iunit, NMRKS(:,I), .false.)
               CALL GETSYM(NMRKS(:,I),NEL,G1,NBASISMAX,ISYM)
               CALL WRITESYM(iunit,ISym%Sym,.TRUE.)
            ENDDO
            CLOSE(iunit)
         endif

!Update: 14.03.2018, K.Ghanem
!FDET has already been assigned in DetPreFreezeInit.
!No idea why it is overwirtten here.
!Therefore, I comment out the following code and hope for the best.
!Instead, I look for FDET in the list of determinants NMRKS and assign the index to IFDET.

!C.. Now generate the fermi determiant
!C.. Work out the fermi det
         !DO I=1,NEL
            !FDET(I)=NMRKS(I,IFDET)
         !ENDDO
         !WRITE(6,*) "Fermi Determinant:",IFDET
         !WRITE(6,*) "Reference determinant to be used for diagonalisation procedure: "
         !call write_det (6, FDET, .true.)

         !if (tDefineDet) then
             !DO I=1,NEL
                 !IF(DefDet(i+NFROZEN)-NFROZEN.ne.FDET(I)) THEN
                     !WRITE(6,"(A)") "*** WARNING - Defined determinant does not match reference determinant in CI matrix ***"
                     !WRITE(6,*) NMRKS(:,IFDET)
                     !WRITE(6,*) DefDet(:)
                     !EXIT
                 !ENDIF
             !ENDDO
         !ENDIF
        IFDET=0
        DO I=1,NDET
            IF(ALL(NMRKS(:,I).EQ.FDET))THEN
                IFDET=I
                Exit
            END IF
        END DO
        IF(IFDET.EQ.0) call stop_all("DetCalcInit","Fermi determinant is not found in NMRKS!")

         WRITE(6,*) ' NUMBER OF SYMMETRY UNIQUE DETS ' , NDET

!C         WRITE(6,*) ' TOTAL NUMBER OF DETS.' , NDETTOT
         IF(NEVAL.EQ.0) THEN
            WRITE(6,*) 'NEVAL=0.  Setting NEVAL=NDET'
            NEVAL=NDET
         ENDIF
         IF(NEVAL.GT.NDET) THEN
            WRITE(6,*) 'NEVAL>NDET.  Setting NEVAL=NDET'
            NEVAL=NDET
         ENDIF

         IF(ABS(DETINV).GT.NDET) THEN
            WRITE(6,*) 'DETINV=',DETINV,'>NDET=',NDET
            WRITE(6,*) 'Setting DETINV to 0'
            DETINV=0
         ENDIF
         CALL neci_flush(6)

!C ==----------------------------------------------------------------==
!C..Set up memory for c's, nrow and the label
         IF(TCALCHMAT) THEN
            WRITE(6,*) "CK Size",NDET*NEVAL*HElement_t_size
            Allocate(CkN(nDet,nEval), stat=ierr)
            LogAlloc(ierr,'CKN',nDet*nEval, HElement_t_sizeB, tagCKN)
            CKN=(0.0_dp)
!C..
            Allocate(Ck(nDet,nEval), stat=ierr)
            LogAlloc(ierr,'CK',nDet*nEval, HElement_t_sizeB, tagCK)
            CK=(0.0_dp)
!C..
            allocate(W(nEval), stat=ierr)
            LogAlloc(ierr, 'W', nEval,8,tagW)
            W=0.0_dp
         ENDIF
!C..
         IF(TREAD) THEN
           CALL READ_PSI(BOX,BOA,COA,NDET,NEVAL,NBASISMAX,NEL,CK,W)
         ENDIF
      ENDIF

      if (tGAS) then
          if (GAS_exc_gen == possible_GAS_exc_gen%GENERAL) then
              call init_GAS()
          else if (GAS_exc_gen == possible_GAS_exc_gen%DISCONNECTED) then
              call init_disconnected_GAS(GAS_specification)
          end if
      end if

!      TMC=TCALCHMAT.AND.(.NOT.TENERGY)

    End Subroutine DetCalcInit

    Subroutine DoDetCalc
      Use global_utilities
      use util_mod, only: get_free_unit
      use Determinants , only : get_helement,FDet, DefDet, tDefineDet
      use SystemData, only : Alat, arr, brr, boa, box, coa, ecore, g1,Beta
      use SystemData, only : nBasis, nBasisMax,nEl,nMsh,LzTot, TSPN,LMS
      use IntegralsData, only: FCK,NMAX, UMat
      Use LoggingData, only: iLogging,tLogDets, tCalcVariationalEnergy
      use SystemData, only  : tCSFOLD
      use Parallel_neci, only : iProcIndex
      use DetBitops, only: DetBitEQ,EncodeBitDet,FindBitExcitLevel
      use bit_rep_data, only: NIfDBO,NIfTot,NIfD
      use legacy_data, only: irat
      use bit_reps, only: decode_bit_det
      use sym_mod
      use HElem
      use MemoryManager, only: TagIntType
      use hist_data, only: tHistSpawn
      use CalcData, only: tFCIDavidson
      !use davidson_neci, only: DavidsonCalcType, DestroyDavidsonCalc
      !use davidson_neci, only: perform_davidson
      !use sparse_arrays, only: calculate_sparse_hamiltonian, deallocate_sparse_ham, &
                               !sparse_ham, hamil_diag, HDiagTag, SparseHamilTags
      !use hamiltonian_linalg, only: sparse_hamil_type
      !use constants, only: size_n_int

      use davidson_neci, only: DavidsonCalcType, DestroyDavidsonCalc
      use davidson_neci, only: davidson_direct_ci_init, davidson_direct_ci_end, perform_davidson
      use hamiltonian_linalg, only: direct_ci_type, tCalcHFIndex
      use FCIMCData, only: davidson_ras, davidson_classes
      use ras, only: generate_entire_ras_space

      real(dp) , ALLOCATABLE :: TKE(:),A(:,:),V(:),AM(:),BM(:),T(:),WT(:),SCR(:),WH(:),WORK2(:),V2(:,:),FCIGS(:)
      HElement_t(dp), ALLOCATABLE :: WORK(:)
      INTEGER , ALLOCATABLE :: INDEX(:),ISCR(:),Temp(:)
      integer(TagIntType) :: TKETag=0,ATag=0,VTag=0,AMTag=0,BMTag=0,TTag=0
      INTEGER(TagIntType) :: WTTag=0,SCRTag=0,ISCRTag=0,INDEXTag=0,WHTag=0,Work2Tag=0,V2Tag=0,WorkTag=0
      integer :: ierr,Lz
      character(25), parameter :: this_routine = 'DoDetCalc'
      real(dp) EXEN,GSEN
      Type(BasisFn) ISym,IHFSYM
      INTEGER GC,I,ICMAX
      INTEGER IN,IND,INDZ
      INTEGER NBLOCK!,OpenOrbs,OpenOrbsSym,Ex(2,NEl)
      INTEGER nKry1
      INTEGER(KIND=n_int) :: ilut(0:NIfTot), ilut_temp(0:NIfTot)
      INTEGER J,JR,iGetExcitLevel_2,ExcitLevel, iunit
      INTEGER LSCR,LISCR,MaxIndex
      LOGICAL tMC!,TestClosedShellDet,Found,tSign
      real(dp) GetHElement, calct, calcmcen, calcdlwdb,norm
      integer:: ic,TempnI(NEl),MomSymDet(NEl),ICSym,ICConnect,PairedUnit,SelfInvUnit
      integer(n_int) :: iLutMomSym(0:NIfTot)
      logical :: tSuccess
      type(DavidsonCalcType) :: davidsonCalc
      integer(n_int), allocatable, dimension(:,:) :: davidson_ilut
      integer, allocatable, dimension(:) :: davidson_parities
      integer :: nI(nel)
      integer :: davidson_size
      !character (len=*), parameter :: t_r = "DoDetCalc"
      !integer(TagIntType) :: IlutTag

      IF(tEnergy) THEN
          WRITE(6,'(1X,A,E19.3)') ' B2LIMIT : ' , B2L
          WRITE(6,*) ' NBLK : ' , NBLK
          WRITE(6,*) ' NKRY : ' , NKRY
          WRITE(6,*) ' NEVAL : ' , NEVAL

          WRITE(6,*) ' NCYCLE : ' , NCYCLE
          WRITE(6,*) ' TENERGY : ' , TENERGY
          WRITE(6,*) ' IOBS : ' , IOBS
          WRITE(6,*) ' JOBS : ' , JOBS
          WRITE(6,*) ' KOBS : ' , KOBS
          WRITE(6,*) ' NMSH : ' , NMSH
          IF(IOBS.GT.NMSH.OR.IOBS.LE.0.OR.JOBS.GT.NMSH.OR.JOBS.LE.0.OR.KOBS.GT.NMSH.OR.KOBS.LE.0) THEN
            call stop_all(this_routine, ' !!! REFERENCE PARTICLE NOT IN BOX !!! ')
          ENDIF
      ENDIF

!C.. now back to the storing H
      IF(TCALCHMAT) THEN
         WRITE(6,*) "Calculating H matrix"
!C..We need to measure HAMIL and LAB first
         ALLOCATE(NROW(NDET),stat=ierr)
         CALL LogMemAlloc('NROW',NDET,4,this_routine,NROWTag,ierr)
         NROW(1:NDET)=0
         ICMAX=1
!Falsify tMC
         TMC=.FALSE.
         allocate(HAMIL(0), LAB(0))
         CALL DETHAM(NDET,NEL,NMRKS,HAMIL,LAB,NROW,.TRUE.,ICMAX,GC,TMC)
         deallocate(HAMIL, LAB)
         WRITE(6,*) ' FINISHED COUNTING '
         WRITE(6,*) "Allocating memory for hamiltonian: ",GC*2
         CALL neci_flush(6)
!C..Now we know size, allocate memory to HAMIL and LAB
         LENHAMIL=GC
         Allocate(Hamil(LenHamil), stat=ierr)
         LogAlloc(ierr, 'HAMIL', LenHamil, HElement_t_sizeB, tagHamil)
         HAMIL=(0.0_dp)
!C..
         ALLOCATE(LAB(LENHAMIL),stat=ierr)
         CALL LogMemAlloc('LAB',LenHamil,4,this_routine,LabTag,ierr)

         LAB(1:LENHAMIL)=0
!C..Now we store HAMIL and LAB
         CALL DETHAM(NDET,NEL,NMRKS,HAMIL,LAB,NROW,.FALSE.,ICMAX,GC,TMC)

         IF(BTEST(ILOGGING,7)) THEN
!C.. we write out H now
            iunit = get_free_unit()
            OPEN(iunit,FILE='HAMIL',STATUS='UNKNOWN')
            J=0
            JR=0
!C            HMAX=-dflmax()
!C            HMIN=dflmax()
            DO I=1,LENHAMIL
               DO WHILE(I.GT.J)
                  JR=JR+1
                  J=J+NROW(JR)
               ENDDO
               WRITE(iunit,"(2I12)",advance='no') JR,LAB(I)
               IF(HElement_t_size.EQ.1) THEN
                  WRITE(iunit,*) HAMIL(I)
               ELSE
                  WRITE(iunit,*) HAMIL(I),ABS(HAMIL(I))
               ENDIF
!C               CALL WRITEDET(14,NMRKS(1,JR),NEL,.FALSE.)
!C               WRITE(14,"(A)",advance='no'),"|"
!C               CALL WRITEDET(14,NMRKS(1,LAB(I)),NEL,.FALSE.)
!C              WRITE(14,"(F27.20)") HAMIL(I)
!C               CALL WRITEDET(14,NMRKS(1,LAB(I)),NEL,.FALSE.)
!C               WRITE(14,"(A)",advance='no'),"|"
!C               CALL WRITEDET(14,NMRKS(1,JR),NEL,.FALSE.)
!C               WRITE(14,"(F27.20)") HAMIL(I)

!C               IF(HAMIL(I).GT.HMAX) HMAX=HAMIL(I)
!C               IF(HAMIL(I).LT.HMIN) HMIN=HAMIL(I)
            ENDDO
            CLOSE(iunit)
         ENDIF
        WRITE(6,*) '<D0|H|D0>=',real(GETHELEMENT(IFDET,IFDET,HAMIL,LAB,NROW,NDET), dp)
        WRITE(6,*) '<D0|T|D0>=',CALCT(NMRKS(1,IFDET),NEL)
        CALL neci_flush(6)
!CC         CALL HAMHIST(HMIN,HMAX,LENHAMIL,NHISTBOXES)
      ENDIF
!C.. We've now finished calculating H if we were going to.
!C.. IF ENERGY CALC (for which we need to have calced H)
!
      IF(TENERGY) THEN
         IF(NBLK.NE.0) THEN
!C..Things needed for Friesner-Pollard diagonalisation
            IF(TMC) call stop_all(this_routine, 'TMC and TENERGY set - Stopping')
            IF(HElement_t_size.NE.1)  call stop_all(this_routine, 'Cannot do Lanczos on Complex orbitals.')
            NKRY1=NKRY+1
            NBLOCK=MIN(NEVAL,NBLK)
            LSCR=MAX(NDET*NEVAL,8*NBLOCK*NKRY)
            LISCR=6*NBLOCK*NKRY
!C..
!            write (6,'(/,/,8X,64(1H*))')
            write (6,'(7X," *",62X,"*")')
          write (6,'(7X," *",19X,A,18X,"*")') ' LANCZOS DIAGONALISATION '
            write (6,'(7X," *",62X,"*")')
!            write (6,'(7X,1X,64(1H*))')
!C..Set up memory for FRSBLKH

            ALLOCATE(A(NEVAL,NEVAL),stat=ierr)
            CALL LogMemAlloc('A',NEVAL**2,8,this_routine,ATag,ierr)
            A=0.0_dp
!C..
!C,, W is now allocated with CK
!C..
            ALLOCATE(V(NDET*NBLOCK*NKRY1),stat=ierr)
            CALL LogMemAlloc('V',NDET*NBLOCK*NKRY1,8,this_routine,VTag,ierr)
            V=0.0_dp
!C..
            ALLOCATE(AM(NBLOCK*NBLOCK*NKRY1),stat=ierr)
            CALL LogMemAlloc('AM',NBLOCK*NBLOCK*NKRY1,8,this_routine,AMTag,ierr)
            AM=0.0_dp
!C..
            ALLOCATE(BM(NBLOCK*NBLOCK*NKRY),stat=ierr)
            CALL LogMemAlloc('BM',NBLOCK*NBLOCK*NKRY,8,this_routine,BMTag,ierr)
            BM=0.0_dp
!C..
            ALLOCATE(T(3*NBLOCK*NKRY*NBLOCK*NKRY),stat=ierr)
            CALL LogMemAlloc('T',3*NBLOCK*NKRY*NBLOCK*NKRY,8,this_routine,TTag,ierr)
            T=0.0_dp
!C..
            ALLOCATE(WT(NBLOCK*NKRY),stat=ierr)
            CALL LogMemAlloc('WT',NBLOCK*NKRY,8,this_routine,WTTag,ierr)
            WT=0.0_dp
!C..
            ALLOCATE(SCR(LScr),stat=ierr)
            CALL LogMemAlloc('SCR',LScr,8,this_routine,SCRTag,ierr)
            SCR=0.0_dp
            ALLOCATE(ISCR(LIScr),stat=ierr)
            CALL LogMemAlloc('IScr',LIScr,4,this_routine,IScrTag,ierr)
            ISCR(1:LISCR)=0
            ALLOCATE(INDEX(NEVAL),stat=ierr)
            CALL LogMemAlloc('INDEX',NEVAL,4,this_routine,INDEXTag,ierr)
            INDEX(1:NEVAL)=0
!C..
            ALLOCATE(WH(NDET),stat=ierr)
            CALL LogMemAlloc('WH',NDET,8,this_routine,WHTag,ierr)
            WH=0.0_dp
            ALLOCATE(WORK2(3*NDET),stat=ierr)
            CALL LogMemAlloc('WORK2',3*NDET,8,this_routine,WORK2Tag,ierr)
            WORK2=0.0_dp
            ALLOCATE(V2(NDET,NEVAL),stat=ierr)
            CALL LogMemAlloc('V2',NDET*NEVAL,8,this_routine,V2Tag,ierr)
            V2=0.0_dp
!C..Lanczos iterative diagonalising routine
            if (t_non_hermitian) then
                call stop_all(this_routine, &
                    "NECI_FRSBLKH not adapted for non-hermitian Hamiltonians")
            end if
            CALL NECI_FRSBLKH(NDET,ICMAX,NEVAL,HAMIL,LAB,CK,CKN,NKRY,NKRY1,NBLOCK,NROW,LSCR,LISCR,A,W,V,AM,BM,T,WT, &
     &  SCR,ISCR,INDEX,NCYCLE,B2L,.true.,.false.,.false.,.true.)

!Multiply all eigenvalues by -1.
            CALL DSCAL(NEVAL,-1.0_dp,W,1)
         ELSE
!C.. We splice in a non-Lanczos diagonalisin routine if NBLOCK=0
            IF(NEVAL.NE.NDET) THEN
               WRITE(6,*) "NEVAL.NE.NDET.",NEVAL,NDET," Cannot exactly diagonalize."
               call stop_all(this_routine, "Cannot exactly diagonalise")
            ENDIF
            WRITE(6,*) "NBLK=0.  Doing exact diagonalization."
            IF(TCALCHMAT) THEN
               ALLOCATE(WORK(4*NDET),stat=ierr)
               CALL LogMemAlloc('WORK',4*NDET,8*HElement_t_size,this_routine,WorkTag,ierr)
               ALLOCATE(WORK2(3*NDET),stat=ierr)
               CALL LogMemAlloc('WORK2',3*NDET,8,this_routine,WORK2Tag,ierr)
               if (t_non_hermitian) then
                   call stop_all(this_routine, &
                       "HDIAG_nec is not setup for non-hermitian Hamiltonians")
               end if
               CALL HDIAG_neci(NDET,HAMIL,LAB,NROW,CK,W,WORK2,WORK,NBLOCKSTARTS,NBLOCKS)
            ENDIF
         ENDIF
!C..
!  Since we no longer use HAMIL or LAB, we deallocate
         if(.not.tCalcVariationalEnergy) then
             LogDealloc(tagHamil)
             Deallocate(Hamil)
             DEALLOCATE(LAB)
             CALL LogMemDealloc(this_routine,LabTag)
         endif
         ALLOCATE(TKE(NEVAL),stat=ierr)
         CALL LogMemAlloc('TKE',NEVAL,8,this_routine,TKETag,ierr)

         EXEN=CALCMCEN(NEVAL,W,BETA)
         WRITE(6,"(A,F19.9)") "EXACT E(BETA)=",EXEN
         GSEN=CALCDLWDB(IFDET,NDET,NEVAL,CK,W,BETA)
         WRITE(6,"(A,F19.9)") "EXACT DLWDB(D0)=",GSEN
         WRITE(6,"(A,F19.9)") "GROUND E=",W(1)
!C.. END ENERGY CALC
!      ENDIF
       ELSEIF(tFCIDavidson)THEN
          if(.not.TSPN .or. LMS.NE.0)then
              call stop_all("DoDetCalc","FCI-Davidson only works for closed shell systems.")
          end if
          davidsonCalc = davidson_direct_ci_init()
          davidson_size = davidsonCalc%super%space_size
          allocate(davidson_ilut(0:NIfTot,davidson_size))
          allocate(davidson_parities(davidson_size))
          call generate_entire_ras_space(davidson_ras, davidson_classes, davidson_size, davidson_ilut, davidson_parities)
          !Find HF index
          !Set this flag, otherwise hfindex will be overwritten
          tCalcHFIndex = .false.
          davidsonCalc%super%hfindex=0
          CALL EncodeBitDet(FDet,iLut(0:NIfDBO))
          do i=1,davidson_size
            if(DetBitEq(davidson_ilut(:,i),ilut))then
                davidsonCalc%super%hfindex=i
                exit
            end if
          end do
          IF(davidsonCalc%super%hfindex.EQ.0) call stop_all("DoDetCalc","Fermi determinant is not found in RAS space!")
          if (t_non_hermitian) then
              call stop_all(this_routine, &
                  "perform_davidson not adapted for non-hermitian Hamiltonians!")
          end if
          call perform_davidson(davidsonCalc, direct_ci_type, .true.)
      ENDIF

      call neci_flush(6)
!C.. If we're calculating rhos (for which we have to have calced H
!No longer used
!      IF(TRHOIJ) THEN
!         IF((.NOT.TENERGY).AND.(.NOT.TREADRHO)) THEN
!            WRITE(6,*) "Calculating approx RHOs"
!            WRITE(6,*) "Using Trotter decomposition? ",TTROT
!            WRITE(6,*) "Order of Taylor: ",ABS(NTAY)
!            CALL CALCAPPROXRHOIJ(BETA,I_P,HAMIL,LAB,NROW,NDET,RHOMIN,RHOMAX,NRHOS,RHOEPS,TTROT,NTAY)
!         ENDIF
!      ENDIF
!C..Free HAMIL AND LAB memory if we no longer need them
!      IF(TCALCHMAT.AND..NOT.(TMONTE.AND.TMC)) THEN
!      ENDIF


!C.. IF we want to compress the found determinants for use later...
      IF(tFindDets) THEN
         IF(tCompressDets) THEN
!Need to find symmetry of the reference determinant, so that we can only look for determinants of the correct symmetry.
            CALL GETSYM(FDET,NEL,G1,NBASISMAX,IHFSYM)
            IF(.not.associated(NMRKS)) THEN
                WRITE(6,*) "NMRKS not allocated"
                CALL neci_flush(6)
                CALL Stop_All("DoDetCalc","NMRKS not allocated so cannot compress dets.")
            ENDIF
!First, we want to count the number of determinants of the correct symmetry...
            Det=0
            norm=0.0_dp
            if(tFCIDavidson)then
                do i=1,davidson_size
                    CALL decode_bit_det(nI, davidson_ilut(:,i))
                    CALL GETSYM(nI,NEL,G1,NBASISMAX,ISYM)
                    !CALL GetLz(nI,NEL,Lz)
                    !IF((ISym%Sym%S.eq.IHFSYM%Sym%S).and.(Lz.eq.LzTot)) THEN
                    IF(ISym%Sym%S.eq.IHFSYM%Sym%S) THEN
                        Det=Det+1
                        norm=norm+(davidsonCalc%davidson_eigenvector(i))**2
                    ENDIF
                enddo
            else
                do i=1,NDET
                    CALL GETSYM(NMRKS(:,i),NEL,G1,NBASISMAX,ISYM)
                    !CALL GetLz(NMRKS(:,i),NEL,Lz)
                    !IF((ISym%Sym%S.eq.IHFSYM%Sym%S).and.(Lz.eq.LzTot)) THEN
                    IF(ISym%Sym%S.eq.IHFSYM%Sym%S) THEN
                        Det=Det+1
                        IF(tEnergy) THEN
                            norm=norm+(REAL(CK(i,1),dp))**2
                        ENDIF
                    ENDIF
                enddo
            end if
            WRITE(6,"(I25,A,I4,A)") Det," determinants of symmetry ",IHFSym%Sym%S," found."
            WRITE(6,*) "Normalization of eigenvector 1 is: ", norm
            CALL neci_flush(6)

            ALLOCATE(FCIDets(0:NIfTot,Det),stat=ierr)
            IF(ierr.ne.0) CALL Stop_All("DetCalc","Cannot allocate memory to hold vector")
            ALLOCATE(FCIDetIndex(0:NEl+1),stat=ierr) !+1 so we can store the end of the array too
            ALLOCATE(Temp(Det),stat=ierr)
            IF(ierr.ne.0) CALL Stop_All("DetCalc","Cannot allocate memory to hold vector")
            IF(tEnergy .or. tFCIDavidson) THEN
                ALLOCATE(FCIGS(Det),stat=ierr)
                IF(ierr.ne.0) CALL Stop_All("DetCalc","Cannot allocate memory to hold vector")
            ENDIF
            if(tCalcVariationalEnergy) then
                !This allows us to resort to get back to the hamiltonian ordering
                allocate(ReIndex(Det),stat=ierr)
                if(ierr.ne.0) CALL Stop_All("DetCalc","Cannot allocate memory to hold vector")
                ReIndex(:)=0
            endif



            Det=0
            FCIDetIndex(:)=0
            if(tFCIDavidson)then
                do i=1,davidson_size
                    CALL decode_bit_det(nI, davidson_ilut(:,i))
                    CALL GETSYM(nI,NEL,G1,NBASISMAX,ISYM)
                    !CALL GetLz(nI,NEL,Lz)
                    !IF((ISym%Sym%S.eq.IHFSYM%Sym%S).and.(Lz.eq.LzTot)) THEN
                    IF(ISym%Sym%S.eq.IHFSYM%Sym%S) THEN
                        Det=Det+1
                        ExcitLevel= FindBitExcitLevel(davidson_ilut(:,i), iLut) !iGetExcitLevel_2(FDet,nI,NEl,NEl)
                        FCIDetIndex(ExcitLevel+1)=FCIDetIndex(ExcitLevel+1)+1
                        Temp(Det)=ExcitLevel    !Temp will now temporarily hold the excitation level of the determinant.
                        FCIDets(:,Det) = davidson_ilut(:, i)
                        FCIGS(Det)=davidson_parities(i)*davidsonCalc%davidson_eigenvector(i)/norm
                    ENDIF
                    if(tCalcVariationalEnergy) ReIndex(i)=i
                enddo
            else
                do i=1,NDet
                    CALL GETSYM(NMRKS(:,i),NEL,G1,NBASISMAX,ISYM)
                    !CALL GetLz(NMRKS(:,i),NEL,Lz)
    !                IF((NMRKS(1,i).eq.28).and.(NMRKS(2,i).eq.29).and.(NMRKS(3,i).eq.30).and.(NMRKS(4,i).eq.31)) THEN
    !                    WRITE(6,*) "Found Det: ",NMRKS(:,i)
    !                    WRITE(6,*) i,iSym%Sym%S,REAL(CK(i),8)
    !                ENDIF
                    !IF((ISym%Sym%S.eq.IHFSYM%Sym%S).and.(Lz.eq.LzTot)) THEN
                    IF(ISym%Sym%S.eq.IHFSYM%Sym%S) THEN
                        Det=Det+1
                        ExcitLevel=iGetExcitLevel_2(FDet,NMRKS(:,i),NEl,NEl)
    ! FCIDetIndex is off by one, for later cumulative indexing
                        FCIDetIndex(ExcitLevel+1)=FCIDetIndex(ExcitLevel+1)+1
                        Temp(Det)=ExcitLevel    !Temp will now temporarily hold the excitation level of the determinant.
                        CALL EncodeBitDet(NMRKS(:,i),FCIDets(0:NIfTot,Det))
                        IF(tEnergy) THEN
                            FCIGS(Det)=REAL(CK(i,1),dp)/norm
                        ENDIF
                    ENDIF
                    if(tCalcVariationalEnergy) ReIndex(i)=i
                enddo
            end if

            IF(iExcitLevel.le.0) THEN
                MaxIndex=NEl
            ELSE
                MaxIndex=MIN(iExcitLevel,NEl)
            ENDIF
!NB FCIDetIndex is off by 1 for later cumulation
            do i=1,MaxIndex
                WRITE(6,*) "Number at excitation level: ",i," is: ",FCIDetIndex(i+1)
            enddo

            ! We now want to sort the determinants according to the
            ! excitation level (stored in Temp)
            IF(.not.tEnergy .and. .not. tFCIDavidson) THEN
                call sort (temp(1:Det), FCIDets(:,1:Det))
            ELSE
                if(tCalcVariationalEnergy) then
                    call sort (temp(1:Det), FCIDets(:,1:Det), FCIGS(1:Det), ReIndex(1:Det))
                else
                    call sort (temp(1:Det), FCIDets(:,1:Det), FCIGS(1:Det))
                endif
!                CALL Stop_All("DetCalc","Cannot do histogramming FCI without JUSTFINDDETS at the
                    !moment (need new sorting - bug ghb24)")
            ENDIF

!Test that HF determinant is the first determinant
            CALL EncodeBitDet(FDet,iLut(0:NIfDBO))
            do i=0,NIfDBO
                IF(iLut(i).ne.FCIDets(i,1)) THEN
                    CALL Stop_All("DetCalc","Problem with ordering the determinants by excitation level")
                ENDIF
            enddo

!Change it so that FCIDetIndex indexes the start of the block of that excitation level.
!            FCIDetIndex(1)=2    !Singles start at index 2
            FCIDetIndex(0)=1
            do i=1,NEl+1
                FCIDetIndex(i)=FCIDetIndex(i-1)+FCIDetIndex(i)
            enddo
            FCIDetIndex(nEl+1)=Det+1
            IF(FCIDetIndex(nEl+1).ne.Det+1) THEN
                WRITE(6,*) "FCIDETIndex:", FCIDetIndex(:)
                CALL Stop_All("DetCalc","Error in the indexing of determinant excitation level")
            ENDIF

!We now need to sort within the excitation level by the "number" of the determinant
            do i=1,MaxIndex
                IF(.not.tEnergy .and. .not. tFCIDavidson) THEN
!                    WRITE(6,*) i,FCIDetIndex(i),FCIDetIndex(i+1)-1
!                    CALL neci_flush(6)
                    call sort (FCIDets(:,FCIDetIndex(i):FCIDetIndex(i+1)-1), &
                               temp(FCIDetIndex(i):FCIDetIndex(i+1)-1))
                ELSE
                    if(tCalcVariationalEnergy) then
                        call sort (FCIDets(:,FCIDetIndex(i):FCIDetIndex(i+1)-1), &
                                   temp(FCIDetIndex(i):FCIDetIndex(i+1)-1), &
                                   FCIGS(FCIDetIndex(i):FCIDetIndex(i+1)-1), &
                                   ReIndex(FCIDetIndex(i):FCIDetIndex(i+1)-1))
                    else
                        call sort (FCIDets(:,FCIDetIndex(i):FCIDetIndex(i+1)-1), &
                                   temp(FCIDetIndex(i):FCIDetIndex(i+1)-1), &
                                   FCIGS(FCIDetIndex(i):FCIDetIndex(i+1)-1))
                     endif
                ENDIF
            enddo

            IF(tEnergy .or. tFCIDavidson) THEN
                IF(tLogDETS.and.iProcIndex.eq.0) THEN
                    iunit = get_free_unit()
                    OPEN(iunit,FILE='SymDETS',STATUS='UNKNOWN')

                    do i=1,Det
                        WRITE(iunit,"(2I17)",advance='no') i,temp(i)
                        do j=0,NIfDBO
                           WRITE(iunit,"(I17)",advance='no') FCIDets(j,i)
                        enddo
                        WRITE(iunit,"(A,G25.16,A)",advance='no') " ",FCIGS(i),"  "
                        Call WriteBitDet(iunit,FCIDets(:,i),.true.)
                   enddo
                   CLOSE(iunit)
                ENDIF
                DEALLOCATE(FCIGS)
            ELSE
                IF(tLogDETS.and.iProcIndex.eq.0) THEN
                    iunit = get_free_unit()
                    OPEN(iunit,FILE='SymDETS',STATUS='UNKNOWN')
                    WRITE(iunit,*) "FCIDETIndex: ",FCIDetIndex(:)
                    WRITE(iunit,*) "***"
                    do i=1,Det
                        WRITE(iunit,"(2I13)",advance='no') i,temp(i)
                        do j=0,NIfDBO
                           WRITE(iunit,"(I13)",advance='no') FCIDets(j,i)
                        enddo
                        WRITE(iunit,"(A)",advance='no') " "
                        Call WriteBitDet(iunit,FCIDets(:,i),.true.)
                    enddo
                    CLOSE(iunit)
                ENDIF
            ENDIF !tEnergy - for dumping compressed ordered GS wavefunction
            DEALLOCATE(Temp)
!            do i=1,Det
!                CALL DecodeBitDet(nK,iLut)
!                CALL GETSYM(nK,NEL,G1,NBASISMAX,ISYM)
!                IF((nK(1).eq.28).and.(nK(2).eq.29).and.(nK(3).eq.30).and.(nK(4).eq.31)) THEN
!                    WRITE(6,*) "Found Det: ",nK(:)
!                    WRITE(6,*) i,iSym%Sym%S,FCIGS(i)
!                ENDIF
!            enddo

!             Det=0
!             maxdet=0
!             do i=1,nel
!                 maxdet=maxdet+2**(nbasis-i)
!             enddo
!             IF(.not.associated(NMRKS)) THEN
!                 WRITE(6,*) "NMRKS not allocated"
!                 CALL neci_flush(6)
!             ENDIF
!             norm=0.0_dp
!             OPEN(17,FILE='SymDETS',STATUS='UNKNOWN')
!
!             do i=1,MAXDET
!                 Bits=0
!                 do j=0,nbasis-1
!                     IF(BTEST(i,j)) THEN
!                         Bits=Bits+1
!                     ENDIF
!                 enddo
!                 IF(Bits.eq.NEl) THEN
!
!                     DO j=1,NDET
!                         CALL EncodeBitDet(NMRKS(:,j),iLut)
!                         IF(iLut(0).eq.i) THEN
!
!                             CALL GETSYM(NMRKS(1,j),NEL,G1,NBASISMAX,ISYM)
!                             IF(ISym%Sym%S.eq.0) THEN
!                                 Det=Det+1
!                                 WRITE(17,"(2I12)",advance='no') Det,iLut(0)
!                                 HEL=GetHElement3(NMRKS(:,j),NMRKS(:,j),0)
!                                 norm=norm+(REAL(CK(j),8))**2
!                                 WRITE(17,"(3G25.16)") REAL(HEL,8),REAL(CK(j),8),norm
!                             ENDIF
!                             EXIT
!                         ENDIF
!                     ENDDO
!                 ENDIF
!
!             ENDDO
!             CLOSE(17)
         ENDIF !tCompressDets
      ENDIF !tFindDets
!C..
      IF(TEnergy) THEN
          IF(.NOT.TCSFOLD) THEN
             CALL CFF_CHCK(NDET,NEVAL,NMRKS,NEL,G1,CK,TKE)
          ELSE
             DO I=1,NEVAL
                TKE(I)=0.0_dp
             ENDDO
          ENDIF
          IF(BTEST(ILOGGING,7)) CALL WRITE_PSI(BOX,BOA,COA,NDET,NEVAL,NBASISMAX,NEL,CK,W)
          IF(BTEST(ILOGGING,8)) CALL WRITE_PSI_COMP(BOX,BOA,COA,NDET,NEVAL,NBASISMAX,NEL,CK,W)
          WRITE(6,*) '       ====================================================== '
          WRITE(6,'(A5,5X,A15,1X,A18,1x,A20)') 'STATE','KINETIC ENERGY', 'COULOMB ENERGY', 'TOTAL ENERGY'
          iunit = get_free_unit()
          OPEN(iunit,FILE='ENERGIES',STATUS='UNKNOWN')
          DO IN=1,NEVAL
             WRITE(6,'(I5,2X,3(F19.11,2x))') IN,TKE(IN),W(IN)-TKE(IN),W(IN)
 !            WRITE(iunit,"(I7)",advance='no') IN
 !            CALL WRITEDET(iunit,NMRKS(1,IN),NEL,.FALSE.)
             WRITE(iunit,"(F19.11)") W(IN)
          ENDDO
          CLOSE(iunit)
          WRITE(6,*)   '       ====================================================== '
      ELSEIF(tFCIDavidson)THEN
          OPEN(iunit,FILE='ENERGIES',STATUS='UNKNOWN')
          WRITE(iunit,"(F19.11)") davidsonCalc%davidson_eigenvalue
          CLOSE(iunit)
!C., END energy calc
      ENDIF

      IF(tFCIDavidson) THEN
          call davidson_direct_ci_end(davidsonCalc)
          call DestroyDavidsonCalc(davidsonCalc)
          DEALLOCATE(davidson_ilut)
          DEALLOCATE(davidson_parities)
      ENDIF
!C.. Jump to here if just read Psi in
      CONTINUE

      !Now deallocate NMRKS if tFindDets and not tEnergy
      if (tFindDets.and.tCompressDets.and.(.not.tEnergy)) then
         DEALLOCATE(NMRKS)
         CALL LogMemDealloc(this_routine,tagNMRKS)
      endif

        End Subroutine DoDetCalc

END MODULE DetCalc

!  Given an exact calculation of eigen-vectors and -values, calculate the expectation value of E(Beta)
      FUNCTION CALCMCEN(NEVAL,W,BETA)
         use constants, only: dp
         IMPLICIT NONE
         INTEGER NEVAL,IK
         real(dp)  W(NEVAL),BETA,DNORM,EN,CALCMCEN
         EN=0.0_dp
         DNORM=0.0_dp
         DO IK=1,NEVAL
            EN=EN+(W(IK))*EXP(-(W(IK)-W(1))*BETA)
            DNORM=DNORM+EXP(-(W(IK)-W(1))*BETA)
         ENDDO
         CALCMCEN=EN/DNORM
         RETURN
      END

!  Given an exact calculation of eigen-vectors and -values, calculate the expectation value of E~(Beta)_I for det I
      FUNCTION CALCDLWDB(I,NDET,NEVAL,CK,W,BETA)
         use constants, only: dp
         IMPLICIT NONE
         INTEGER NDET,NEVAL,IK,I
         HElement_t(dp) CK(NDET,NEVAL)
         real(dp)  W(NEVAL),BETA,DNORM,EN,CALCDLWDB
         EN=0.0_dp
         DNORM=0.0_dp
         DO IK=1,NEVAL
            EN=EN+abs(CK(I,IK))**2*(W(IK))*EXP(-(W(IK)-W(1))*BETA)
            DNORM=DNORM+abs(CK(I,IK))**2*EXP(-(W(IK)-W(1))*BETA)
         ENDDO
         CALCDLWDB=EN/DNORM
         RETURN
      END

      SUBROUTINE CFF_CHCK(NDET,NEVAL,NM,NEL,G1,CG,TKE)
      use constants, only: dp
      use util_mod, only: get_free_unit
      USE OneEInts, only : GetTMATEl
      use SystemData, only: BasisFN
      use HElem
      IMPLICIT NONE
      INTEGER NEL,NM(NEL,*),NDET,NEVAL, iunit
      HElement_t(dp) CG(NDET,NEVAL)
      real(dp) TKE(NEVAL)
      TYPE(BASISFN) G1(*)
      real(dp) PI,S,SUM1
      real(dp) AUX
      INTEGER I,J,IN,IEL,L
!..Calculate the expectation value of the kinetic energy
!..<Psi|T|Psi>
      PI=ACOS(-1.0_dp)
      DO IN=1,NEVAL
        TKE(IN)=0.0_dp
        DO I=1,NDET
          SUM1=0.0_dp
          DO J=1,NEL
            AUX=GetTMATEl(NM(J,I),NM(J,I))
!((ALAT(1)**2)*((G1(1,NM(J,I))**2)/(ALAT(1)**2)+
!     &     (G1(2,NM(J,I))**2)/(ALAT(2)**2)+
!     &     (G1(3,NM(J,I))**2)/(ALAT(3)**2)))
            SUM1=SUM1+(AUX)
          ENDDO
!..Cube multiplier
!          CST=PI*PI/(2.0_dp*ALAT(1)*ALAT(1))
!.. Deal with the UEG
!          IF(NBASISMAX(1,1).LE.0) CST=CST*4.0_dp
!          SUM1=CST*SUM1
          TKE(IN)=TKE(IN)+SUM1*abs(CG(I,IN))**2
        ENDDO
      ENDDO
! ==--------------------------------------------------------------==
      IF(.FALSE.) THEN
!      IF(BTEST(ILOGGING,7)) THEN
      iunit = get_free_unit()
      OPEN(iunit,FILE='PSI',STATUS='UNKNOWN')
      DO J=1,NEVAL
        IF(J.EQ.1) THEN
          WRITE(iunit,*) ' GROUND STATE COEFFICIENTS  '
        ELSE
          WRITE(iunit,*) ' COEFFICIENTS FOR EXCITED STATE NUMBER : ' , J
        ENDIF
        S=0.0_dp
        DO I=1,NDET
         IF(abs(CG(I,J)).gt.1.0e-15_dp) THEN
            DO IEL=1,NEL
               WRITE(iunit,"(I3,I3,2I3,2X)",advance='no') (G1(NM(1,IEL))%K(L),L=1,5)
            ENDDO
            IF(HElement_t_size.EQ.1) THEN
               WRITE(iunit,"(F19.9,1X,I7)") CG(I,J),I
            ELSE
               WRITE(iunit,"(F19.9,1X,I7)") CG(I,J),I
            ENDIF
         ENDIF
         S=S+abs(CG(I,J))**2
        ENDDO
        WRITE(iunit,'(A,F19.5)') ' SQUARE OF COEFFICIENTS : ' , S
        WRITE(iunit,*)
      ENDDO
      CLOSE(iunit)
      ENDIF
      RETURN
      END
