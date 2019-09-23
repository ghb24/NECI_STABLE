module read_fci

    character(len=1024) :: FCIDUMP_name

contains

    SUBROUTINE INITFROMFCID(NEL,NBASISMAX,LEN,LMS,TBIN)
         use SystemData , only : tNoSymGenRandExcits,lNoSymmetry,tROHF,tHub,tUEG
         use SystemData , only : tStoreSpinOrbs,tKPntSym,tRotatedOrbsReal,tFixLz,tUHF
         use SystemData , only : tMolpro,tReadFreeFormat,tReltvy, nclosedOrbs, nOccOrbs
         use SystemData , only : nIrreps
         use SymData, only: nProp, PropBitLen, TwoCycleSymGens
         use Parallel_neci
         use util_mod, only: get_free_unit, near_zero
         IMPLICIT NONE
         logical, intent(in) :: tbin
         integer, intent(out) :: nBasisMax(5,*),LEN,LMS
         integer, intent(in) :: NEL
         integer SYMLZ(1000)
         integer OCC(nIrreps), CLOSED(nIrreps), FROZEN(nIrreps)
         integer(int64) :: ORBSYM(1000)
         INTEGER NORB,NELEC,MS2,ISYM,i,SYML(1000), iunit,iuhf
         LOGICAL exists
         logical :: uhf,trel,tDetectSym

#ifdef _MOLCAS_
         logical :: tExists     !test for existence of input file.
         integer :: isfreeunit  !function returning integer for free unit
#endif

         CHARACTER(len=3) :: fmat
         NAMELIST /FCI/ NORB,NELEC,MS2,ORBSYM,OCC,CLOSED,FROZEN,&
              ISYM,IUHF,UHF,TREL,SYML,SYMLZ,PROPBITLEN,NPROP
         UHF=.FALSE.
         fmat='NO'
         PROPBITLEN = 0
         NPROP = 0
         IUHF = 0
         TREL = .false.
         SYMLZ(:) = 0
         OCC = -1
         CLOSED = -1
         FROZEN = -1
         ! [W.D. 15.5.2017:]
         ! with the new relativistic calculations, withoug a ms value in the
         ! FCIDUMP, we have to set some more defaults..
         MS2 = 0
         IF(iProcIndex.eq.0) THEN
#ifdef _MOLCAS_
           call f_Inquire('FCIDMP',tExists)
           if(tExists) then
             iUnit=17
             iUnit=IsFreeUnit(iUnit)
             Call Molcas_Open(iunit,'FCIDMP')
             Rewind(iunit)
             READ(iunit,FCI)
!             write(6,*) 'FCI NAMELIST print 1'
!             WRITE(6,FCI)
           else
              call Stop_All('InitFromFCID','FCIDUMP file does not exist')
           end if
#else
             iunit = get_free_unit()
             IF(TBIN) THEN
                INQUIRE(FILE='FCISYM',EXIST=exists)
                IF(.not.exists) THEN
                    CALL Stop_All('InitFromFCID','FCISYM file does not exist')
                ENDIF
                INQUIRE(FILE=FCIDUMP_name,EXIST=exists,FORMATTED=fmat)
                IF(.not.exists) THEN
                    CALL Stop_All('INITFROMFCID','FCIDUMP file does not exist')
                ENDIF
                OPEN(iunit,FILE='FCISYM',STATUS='OLD',FORM='FORMATTED')
                READ(iunit,FCI)
             ELSE
                INQUIRE(FILE=FCIDUMP_name,EXIST=exists,UNFORMATTED=fmat)
                IF(.not.exists) THEN
                    CALL Stop_All('InitFromFCID','FCIDUMP file does not exist')
                ENDIF
                OPEN(iunit,FILE=FCIDUMP_name,STATUS='OLD',FORM='FORMATTED')
                READ(iunit,FCI)
             ENDIF
#endif
             CLOSE(iunit)
         ENDIF

!Now broadcast these values to the other processors
         CALL MPIBCast(NORB,1)
         CALL MPIBCast(NELEC,1)
         CALL MPIBCast(MS2,1)
         CALL MPIBCast(ORBSYM,1000)
         CALL MPIBCast(SYML,1000)
         CALL MPIBCast(SYMLZ,1000)
         CALL MPIBCast(ISYM,1)
         CALL MPIBCast(UHF)
         CALL MPIBCast(tRel)
         CALL MPIBCast(PROPBITLEN,1)
         CALL MPIBCast(NPROP,3)
         call MPIBCast(OCC,8)
         call MPIBCast(CLOSED,nIrreps)
         call MPIBCast(FROZEN,nIrreps)
         ! If PropBitLen has been set then assume we're not using an Abelian
         ! symmetry group which has two cycle generators (ie the group has
         ! complex representations).
         TwoCycleSymGens = PropBitLen == 0
         tReltvy = tRel
         nOccOrbs = OCC
         nClosedOrbs = CLOSED

         IF(.not.TwoCycleSymGens.and.((NPROP(1)+NPROP(2)+NPROP(3)).gt.3)) THEN
             !We are using abelian k-point symmetry. Turn it on.
             tKPntSym=.true.
             WRITE(6,"(A,I5,A)") "Using abelian k-point symmetry - ",NPROP(1)*NPROP(2)*NPROP(3)," kpoints found."
         ELSE
             tKPntSym=.false.
         ENDIF

#ifndef __CMPLX
         if(PropBitLen.ne.0) then
             !We have compiled the code in real, but we are looking at a complex FCIDUMP, potentially
             !even with multiple k-points. However, these orbitals must be real, and ensure that the
             !integrals are read in as real.
             write(6,"(A)") "Neci compiled in real mode, but kpoints detected."
             write(6,"(A)") "This will be ok if kpoints are all at Gamma or BZ boundary and are correctly rotated."
             tRotatedOrbsReal=.true.
         elseif(tRel) then
             write(6,"(A)") "Relativistic integrals, but neci compiled in real mode."
         endif
#endif

         if(.not.tMolpro) then
             IF(tROHF.and.(.not.UHF)) THEN
                 CALL Stop_All("INITFROMFCID","ROHF specified, but FCIDUMP is not in a high-spin format.")
             ENDIF
         endif

         tDetectSym=.true.
         DO i=1,NORB
             IF(ORBSYM(i).eq.0.and.TwoCycleSymGens) THEN
                 WRITE(6,*) "** WARNING **"
                 WRITE(6,*) "** Unconverged symmetry of orbitals **"
                 WRITE(6,*) "** Turning point group symmetry off for rest of run **"
                 if(.not.(tFixLz.or.tKPntSym.or.tUEG.or.tHub)) then
                     WRITE(6,*) "** No symmetry at all will be used in excitation generators **"
                     tNoSymGenRandExcits=.true.
                 endif
                 lNoSymmetry=.true.
                 EXIT
             ENDIF
         ENDDO
         if(tDetectSym) then
             !Check that there is more than one irrep
             tDetectSym=.false.
             do i = 1,Norb
                 if(OrbSym(i).ne.1) then
                     tDetectSym=.true.
                     exit
                 endif
             enddo
         endif
         if(.not.tDetectSym) then
             write(6,*) "Only one irrep found. Turning off symmetry for rest of calculation."
             tNoSymGenRandExcits = .true.
             lNoSymmetry = .true.
         endif

         IF(NELEC.NE.NEL) THEN
             WRITE(6,*)                                                 &
     &      '*** WARNING: NEL in FCIDUMP differs from input file ***'
            WRITE(6,*) ' NUMBER OF ELECTRONS : ' , NEL
         ENDIF
!         NEL=NELEC
         IF(LMS.NE.MS2) THEN
            WRITE(6,*)                                                  &
     &      '*** WARNING: LMS in FCIDUMP differs from input file ***'
            LMS=MS2
            WRITE(6,*) ' BASIS MS : ' , LMS
         ENDIF
         if(tMolpro) then
             !Molpros FCIDUMPs always indicate the number of spatial orbitals.
             !Need to x2 to get spin orbitals
             LEN=2*NORB
         else
             IF(UHF.or.tRel) then
                LEN=NORB
             ELSE
                LEN=2*NORB
             ENDIF
         endif
         NBASISMAX(1:5,1:3)=0
!         NBASISMAX(1,1)=1
!         NBASISMAX(1,2)=NORB
         NBASISMAX(1,1)=0
         NBASISMAX(1,2)=0
         NBASISMAX(1,3)=2
         NBASISMAX(4,1)=-1
         NBASISMAX(4,2)=1

!the indicator of UHF (which becomes ISpinSkip or ISS later
         if(tMolpro) then
             if(tUHF) then
                NBASISMAX(2,3)=1
                tStoreSpinOrbs=.true.   !indicate that we are storing the orbitals in umat as spin-orbitals
            else
                NBASISMAX(2,3)=2
                tStoreSpinOrbs=.false.   !indicate that we are storing the orbitals in umat as spatial-orbitals
            endif
         else
             IF((UHF.and.(.not.tROHF)).or.tRel) then
                NBASISMAX(2,3)=1
                tStoreSpinOrbs=.true.   !indicate that we are storing the orbitals in umat as spin-orbitals
             else
                NBASISMAX(2,3)=2
                tStoreSpinOrbs=.false.   !indicate that we are storing the orbitals in umat as spatial-orbitals
             endif
         endif

! Show that there's no momentum conservation
         NBASISMAX(3,3)=1
         RETURN
      END SUBROUTINE INITFROMFCID


      SUBROUTINE GETFCIBASIS(NBASISMAX,ARR,BRR,G1,LEN,TBIN)
         use SystemData, only: BasisFN,BasisFNSize,Symmetry,NullBasisFn,tMolpro,tUHF
         use SystemData, only: tCacheFCIDUMPInts,tROHF,tFixLz,iMaxLz,tRotatedOrbsReal
         use SystemData, only: tReadFreeFormat,SYMMAX,tReltvy, irrepOrbOffset, nIrreps
         use UMatCache, only: nSlotsInit,CalcNSlotsInit
         use UMatCache, only: GetCacheIndexStates,GTID
         use SymData, only: nProp, PropBitLen, TwoCycleSymGens
         use Parallel_neci
         use constants, only: dp,sizeof_int
         use util_mod, only: get_free_unit
         IMPLICIT NONE
         integer, intent(in) :: LEN
         integer, intent(inout) :: nBasisMax(5,*)
         integer, intent(out) :: BRR(LEN)
         real(dp), intent(out) :: ARR(LEN,2)
         type(BasisFN), intent(out) :: G1(LEN)
         HElement_t(dp) Z
         COMPLEX(dp) :: CompInt
         integer(int64) IND,MASK
         INTEGER I,J,K,L,I1
         INTEGER ISYMNUM,ISNMAX,SYMLZ(1000), iunit
         INTEGER NORB,NELEC,MS2,ISYM,ISPINS,ISPN,SYML(1000)
         integer OCC(nIrreps), CLOSED(nIrreps), FROZEN(nIrreps)
         integer(int64) ORBSYM(1000)
         INTEGER nPairs,iErr,MaxnSlot,MaxIndex,IUHF
         INTEGER , ALLOCATABLE :: MaxSlots(:)
         character(len=*), parameter :: t_r='GETFCIBASIS'
         LOGICAL TBIN
         logical :: uhf,tRel
         integer :: orbsPerIrrep(nIrreps)
         NAMELIST /FCI/ NORB,NELEC,MS2,ORBSYM,OCC,CLOSED,FROZEN,&
              ISYM,IUHF,UHF,TREL,SYML,SYMLZ,PROPBITLEN,NPROP

#ifdef _MOLCAS_
         logical :: tExists     !test for existence of input file.
         integer :: isfreeunit  !function returning integer for free unit
#endif
         iunit = 0
         UHF=.FALSE.
         PROPBITLEN = 0
         NPROP = 0
         IUHF = 0
         TREL = .false.
         SYMLZ(:) = 0
         IF(iProcIndex.eq.0) THEN
#ifdef _MOLCAS_
           call f_Inquire('FCIDMP',tExists)
           if(tExists) then
             iUnit=17
             iUnit=IsFreeUnit(iUnit)
             Call Molcas_Open(iunit,'FCIDMP')
             Rewind(iunit)
             READ(iunit,FCI)
!            write(6,*) 'FCI NAMELIST print 2'
!            WRITE(6,FCI)
           else
              call Stop_All('InitFromFCID','FCIDUMP file does not exist')
           end if
#else
             iunit = get_free_unit()
             IF(TBIN) THEN
                OPEN(iunit,FILE='FCISYM',STATUS='OLD',FORM='FORMATTED')
                READ(iunit,FCI)
                CLOSE(iunit)
                OPEN(iunit,FILE=FCIDUMP_name,STATUS='OLD',FORM='UNFORMATTED')
             ELSE
                OPEN(iunit,FILE=FCIDUMP_name,STATUS='OLD',FORM='FORMATTED')
                READ(iunit,FCI)
             ENDIF
#endif
         ENDIF

!Now broadcast these values to the other processors (the values are only read in on root)
         CALL MPIBCast(NORB,1)
         CALL MPIBCast(NELEC,1)
         CALL MPIBCast(MS2,1)
         CALL MPIBCast(ORBSYM,1000)
         CALL MPIBCast(SYML,1000)
         CALL MPIBCast(SYMLZ,1000)
         CALL MPIBCast(ISYM,1)
         CALL MPIBCast(IUHF,1)
         CALL MPIBCast(UHF)
         CALL MPIBCast(TREL)
         CALL MPIBCast(PROPBITLEN,1)
         CALL MPIBCast(NPROP,3)
         ! If PropBitLen has been set then assume we're not using an Abelian
         ! symmetry group which has two cycle generators (ie the group has
         ! complex representations).
         TwoCycleSymGens = PropBitLen == 0
         tReltvy = tRel

         ISYMNUM=0
         ISNMAX=0
         ISPINS=2

         IF((UHF.and.(.not.tROHF)).or.tReltvy) ISPINS=1

         IF(tROHF) THEN
             if(.not.tMolpro) then
                 WRITE(6,*) "Reading in Spin-orbital FCIDUMP but storing as spatial orbitals."
!             WRITE(6,*) "*** Warning - Fock energy of orbitals will be "&
!     &                //"incorrect ***"
!We are reading in the symmetry of the orbitals in spin-orbitals - we need to change this to spatial orbitals
                 DO i=1,NORB-1,2
                     IF(ORBSYM(i).ne.ORBSYM(i+1)) THEN
                         CALL Stop_All("GetFCIBASIS","Spin-orbitals are not ordered in symmetry pairs")
                     ENDIF
                     ORBSYM((i+1)/2)=ORBSYM(i)
                     SYML((i+1)/2)=SYML(i)
                     SYMLZ((i+1)/2)=SYMLZ(i)
                 ENDDO
!NORB is the length in spatial orbitals, LEN is spin orbitals (both spin for UHF)
                 NORB=NORB/2
                 DO i=NORB+1,1000
                     ORBSYM(i)=0
                     SYML(i)=0
                     SYMLZ(i)=0
                 ENDDO
             else
                 !Molpro goes here. tROHF is handled exactly the same way as RHF.
                 !Therefore, the arrays should be the correct length anyway
                 write(6,*) "Reading in ROHF FCIDUMP, and storing as spatial orbitals."
             endif
         ENDIF

         ! Count the number of orbs per irrep
         if(any(ORBSYM(1:NORB).eq.0)) then
            write(6,*) "WARNING: Invalid ORBSYM in FCIDUMP, are you sure you know what you are doing?"
         else
            orbsPerIrrep = 0
            do i = 1, NORB
               orbsPerIrrep(ORBSYM(i)) = orbsPerIrrep(ORBSYM(i)) + ISPINS
            end do
            irrepOrbOffset(1) = 0
            do i = 2, nIrreps
               irrepOrbOffset(i) = irrepOrbOffset(i-1) + orbsPerIrrep(i-1)
            end do
         endif

         !Below is a mistake - since ISPINS is always two, it will be doubled up for us automatically for G1
!         if(tMolpro.and.tUHF) then
!             !Double up the arrays to take into account both spins
!             do i=NORB*2,1,-1
!                 if(mod(i,2).eq.0) then
!                     ORBSYM(i)=ORBSYM(i/2)
!                     SYML(i)=SYML(i/2)
!                     SYMLZ(i)=SYMLZ(i/2)
!                 else
!                     ORBSYM(i)=ORBSYM((i+1)/2)
!                     SYML(i)=SYML((i+1)/2)
!                     SYMLZ(i)=SYMLZ((i+1)/2)
!                 endif
!             enddo
!         endif

         !For Molpro, ISPINS should always be 2, and therefore NORB is spatial, and len is spin orbtials
         IF(LEN.NE.ISPINS*NORB) call stop_all(t_r, 'LEN .NE. NORB in GETFCIBASIS')
         G1(1:LEN)=NullBasisFn
         ARR=0.0_dp

!If we are reading in and cacheing the FCIDUMP integrals, we need to know the maximum number j,l pairs of
!integrals for a given i,k pair. (Or in chemical notation, the maximum number of k,l pairs for a given i,j pair.)
         IF(tCacheFCIDUMPInts) THEN
             WRITE(6,*) "Calculating number of slots needed for "       &
     &           //"integral Cache..."
             IF(tROHF) THEN
                 WRITE(6,*) "This will be inefficient, since multiple " &
     &           //"versions of the same integral are present since we "&
     &           //"are converting to spatial orbitals, so nSlots will "&
     &           //"be too large."
                 CALL Stop_All("readint","ROHF doesn't seem to be "     &
     &            //"working with Caching. ghb can fix if needed...")
             ENDIF
             nPairs=NORB*(NORB+1)/2
!             WRITE(6,*) "NPAIRS: ",NORB,NPAIRS
             CALL neci_flush(6)
             ALLOCATE(MaxSlots(nPairs),stat=ierr)
             MaxSlots(:)=0
         ENDIF


         IF(iProcIndex.eq.0.and.(.not.tMolpro)) THEN
!Just read in integrals on head node. This is only trying to read in fock energies, which aren't written out by molpro anyway

             IF(TBIN) THEN
                if(tMolpro) call stop_all(t_r,'UHF Bin read not functional through molpro')
                IF(UHF.or.tROHF) call stop_all(t_r, 'UHF Bin read not functional')
                MASK=(2**16)-1

                !IND contains all the indices in an integer(int64) - use mask of 16bit to extract them
2               READ(iunit,END=99) Z,IND
                L=int(iand(IND,MASK),sizeof_int)
                IND=Ishft(IND,-16)
                K=int(iand(IND,MASK),sizeof_int)
                IND=Ishft(IND,-16)
                J=int(iand(IND,MASK),sizeof_int)
                IND=Ishft(IND,-16)
                I=int(iand(IND,MASK),sizeof_int)

!                I=Index(I)
!                J=Index(J)
!                K=Index(K)
!                L=Index(L)

!.. Each orbital in the file corresponds to alpha and beta spinorbitals
             !Fill ARR with the energy levels
                IF(I.NE.0.AND.K.EQ.0.AND.I.EQ.J) THEN
                    IF(I.GT.1) THEN
                        IF(ORBSYM(I).NE.ORBSYM(I-1)) THEN
                            IF(ISYMNUM.GT.ISNMAX) ISNMAX=ISYMNUM
                            ISYMNUM=0
                        ENDIF
                    ENDIF
                    ISYMNUM=ISYMNUM+1
                    ARR(2*I-1,1)=real(Z,dp)
                    ARR(2*I,1)=real(Z,dp)
                ELSEIF(I.NE.0.AND.K.EQ.0.AND.J.EQ.0) THEN
                    ARR(2*I-1,1)=real(Z,dp)
                    ARR(2*I,1)=real(Z,dp)
                ENDIF
!.. At the moment we're ignoring the core energy
                IF(I.NE.0) GOTO 2

             ELSE   !Reading in formatted FCIDUMP file

                 ! Can't use * as need to be backward compatible with existing
                 ! FCIDUMP files, some of which have more than 100 basis
                 ! functions and so the integer labels run into each other.
                 ! This means it won't work with more than 999 basis
                 ! functions...
#ifdef __CMPLX
1               READ(iunit,*,END=99) Z,I,J,K,L
#else
1               CONTINUE
                !It is possible that the FCIDUMP can be written out in complex notation, but still only
                !have real orbitals. This occurs with solid systems, where all kpoints are at the
                !gamma point or BZ boundary and have been appropriately rotated. In this case, all imaginary
                !components should be zero to numerical precision.
                if(tRotatedOrbsReal) then
                    !Still need to read in integrals as complex numbers - this FCIDUMP will be from VASP.
                    READ(iunit,*,END=99) CompInt,I,J,K,L
                    Z=real(CompInt,dp)
                    if(abs(aimag(CompInt)).gt.1.0e-7_dp) then
                        write(6,*) "Using the *real* neci compiled code, &
                                   &however non-zero imaginary parts of &
                                   &integrals found"
                        write(6,*) "If all k-points are at Gamma or BZ &
                                   &boundary, orbitals should be able to be &
                                   &rotated to be real"
                        write(6,*) "Check this is the case, or rerun in &
                                   &complex mode to handle complex integrals."
                        call stop_all("GETFCIBASIS", "Real orbitals indicated,&
                                      & but imaginary part of integrals larger&
                                      & than 1.0e-7_dp")
                    endif
                else
                    if(tMolpro.or.tReadFreeFormat) then
                        !If calling from within molpro, integrals are written out to greater precision
                        !Here is where we read integrals in Molcas/NECI interface
                        read(iunit,*,END=99) Z,I,J,K,L
                    else
                        READ(iunit,'(1X,G20.12,4I3)',END=99) Z,I,J,K,L
                    endif
                endif
#endif

                IF(tROHF.and.(.not.tMolpro)) THEN
!The FCIDUMP file is in spin-orbitals - we need to transfer them to spatial orbitals (unless from molpro, where already spatial).
                    IF(I.ne.0) THEN
                        I1 = GTID(I)
                    ELSE
                        I1=0
                    ENDIF
                    IF(J.ne.0) THEN
                        J = GTID(J)
                    ENDIF
                    IF(K.ne.0) THEN
                        K = GTID(K)
                    ENDIF
                    IF(L.ne.0) THEN
                        L = GTID(L)
                    ENDIF
                ELSE
                    I1=I    !Create new I index, since ROHF wants both spin and spatial indicies to get the fock energies right.
                ENDIF


                IF(tCacheFCIDUMPInts) THEN
                    CALL CalcNSlotsInit(I1,J,K,L,Z,nPairs,MaxSlots)
                ENDIF

!.. Each orbital in the file corresponds to alpha and beta spinorbitals
             !Fill ARR with the energy levels
                IF(I1.NE.0.AND.K.EQ.0.AND.I1.EQ.J) THEN
                    IF(I1.GT.1) THEN
                        IF(ORBSYM(I1).NE.ORBSYM(I1-1)) THEN
                            IF(ISYMNUM.GT.ISNMAX) ISNMAX=ISYMNUM
                            ISYMNUM=0
                        ENDIF
                    ENDIF
                    ISYMNUM=ISYMNUM+1

!This fills the single particle energy array (ARR) with the diagonal one-electron integrals, so that if
!there are no fock energies printed out in the FCIDUMP, then we can still order the orbitals in some way.
!If the fock energies are printed out before the one-electron integrals, this will cause problems.
!These integrals must be real.
                    IF(ISPINS.eq.1) THEN
                        ARR(I1,1)=real(Z,dp)
                    ELSEIF(tROHF) THEN
                        ARR(I,1)=real(Z,dp)
                    ELSE
                        ARR(2*I1,1)=real(Z,dp)
                        ARR(2*I1-1,1)=real(Z,dp)
                    ENDIF

                ELSEIF(I1.NE.0.AND.K.EQ.0.AND.J.EQ.0) THEN
!                    WRITE(6,*) I
                    IF(ISPINS.eq.1) THEN
                        ARR(I1,1)=real(Z,dp)
                    ELSEIF(tROHF) THEN
                        ARR(I,1)=real(Z,dp)
                    ELSE
                        ARR(2*I1,1)=real(Z,dp)
                        ARR(2*I1-1,1)=real(Z,dp)
                    ENDIF

!                    DO ISPN=1,ISPINS
!                        ARR(ISPINS*I-ISPN+1,1)=real(Z,dp)
!                    ENDDO

                ENDIF
!.. At the moment we're ignoring the core energy
                IF(I1.NE.0) GOTO 1
             ENDIF
99           CONTINUE
             CLOSE(iunit)

         ENDIF

         if(tMolpro.and.(iProcIndex.eq.0)) close(iunit)

!We now need to broadcast all the information we've just read in...
         IF(tCacheFCIDUMPInts) THEN
             CALL MPIBCast(MaxSlots,nPairs)
         ENDIF
         CALL MPIBCast(ISNMAX,1)
         CALL MPIBCast(ISYMNUM,1)
         CALL MPIBCast(Arr,LEN*2)

         SYMMAX=1
         iMaxLz=0
         DO I=1,NORB
            DO ISPN=1,ISPINS
                BRR(ISPINS*I-ISPN+1)=ISPINS*I-ISPN+1
                if (TwoCycleSymGens) then
                    G1(ISPINS*I-ISPN+1)%Sym%s=ORBSYM(I)-1
                else
                    G1(ISPINS*I-ISPN+1)%Sym%s=ORBSYM(I)
                end if
                IF(tFixLz) THEN
                    G1(ISPINS*I-ISPN+1)%Ml=SYMLZ(I)
                ELSE
                    G1(ISPINS*I-ISPN+1)%Ml=0
                ENDIF
!.. set momentum to 0
                G1(ISPINS*I-ISPN+1)%k(1)=0
                G1(ISPINS*I-ISPN+1)%k(2)=0
                G1(ISPINS*I-ISPN+1)%k(3)=0
                G1(ISPINS*I-ISPN+1)%Ms=-MOD(ISPINS*I-ISPN+1,2)*2+1
                IF(SYMMAX.lt.ORBSYM(I)) SYMMAX=int(ORBSYM(I),sizeof_int)
                IF(abs(SYMLZ(I)).gt.iMaxLz) iMaxLz=abs(SYMLZ(I))
            ENDDO
         ENDDO
         IF(.not.tFixLz) iMaxLz=0
         if (.not.TwoCycleSymGens) then
             SYMMAX = ISYM
         else
             ! We use bit strings to store symmetry information.
             ! SYMMAX needs to be the smallest power of 2 greater or equal to
             ! the actual number of symmetry representations spanned by the basis.
             SYMMAX = 2**ceiling(log(real(SYMMAX,dp))/log(2.0_dp))
         endif
         IF(tFixLz) WRITE(6,"(A,I3)") "Maximum Lz orbital: ",iMaxLz
         WRITE(6,"(A,I3)") "  Maximum number of symmetries: ",SYMMAX
         NBASISMAX(1,1)=0
         NBASISMAX(1,2)=0
         NBASISMAX(5,1)=0
         NBASISMAX(5,2)=SYMMAX-1
         NBASISMAX(2,1)=0
         NBASISMAX(2,2)=0
         IF(tCacheFCIDUMPInts) THEN
!Calculate the minimum required value for nSlotsInit
             MaxnSlot=0
             do i=1,nPairs
                 IF(MaxSlots(i).gt.MaxnSlot) THEN
                     MaxnSlot=MaxSlots(i)
                     MaxIndex=i
                 ENDIF
             enddo

             DEALLOCATE(MaxSlots)
             nSlotsInit=MaxnSlot+1
             WRITE(6,*) "Maximum number of slots needed in cache is: ",nSlotsInit
             WRITE(6,*) "Index corresponding to this maximum number:",MaxIndex
             CALL GetCacheIndexStates(MaxIndex,i,k)
             WRITE(6,"(A,2I6)") "This corresponds to an (i,k) pair of ",i,k
         ENDIF
!         WRITE(6,*) Arr(:,1)
         RETURN
      END SUBROUTINE GETFCIBASIS

!tReadFreezeInts only matters when we are cacheing the FCIDUMP file.
!It is set if we want to cache the integrals to enable the freezing routine to take place, i.e. the <ij|kj> integrals.
!The UMAT2D integrals will also be read in in this case.
!If tReadFreezeInts is false, then if we are cacheing the FCIDUMP file, then we will read and cache all the integrals.
      SUBROUTINE READFCIINT(UMAT,umat_win,NBASIS,ECORE,tReadFreezeInts)
         use constants, only: dp,sizeof_int
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB,NEl
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB,tMolpro
         use SystemData, only: UMatEps,tCacheFCIDUMPInts,tUHF
         use SystemData, only: tRIIntegrals,nBasisMax,tROHF,tRotatedOrbsReal
         use SystemData, only: tReadFreeFormat, G1, tFixLz, tReltvy, nIrreps
         USE UMatCache, only: UMatInd,UMatConj,UMAT2D,TUMAT2D,nPairs,CacheFCIDUMP
         USE UMatCache, only: FillUpCache,GTID,nStates,nSlots,nTypes
         USE UMatCache, only: UMatCacheData,UMatLabels,GetUMatSize
         use OneEInts, only: TMatind,TMat2D,TMATSYM
         use OneEInts, only: CalcTMatSize
         use Parallel_neci
         use shared_memory_mpi
         use SymData, only: nProp, PropBitLen, TwoCycleSymGens
         use util_mod, only: get_free_unit, near_zero
         IMPLICIT NONE
         integer, intent(in) :: NBASIS
         logical, intent(in) :: tReadFreezeInts
         real(dp), intent(out) :: ECORE
         HElement_t(dp), intent(out) :: UMAT(:)
         integer(MPIArg) :: umat_win
         HElement_t(dp) Z
         COMPLEX(dp) :: CompInt
         INTEGER(int64) :: ZeroedInt,NonZeroInt, LzDisallowed
         INTEGER I,J,K,L,X,Y,iSpinType, iunit
         INTEGER NORB,NELEC,MS2,ISYM,SYML(1000)
         integer(int64) ORBSYM(1000)
         LOGICAL LWRITE
         logical :: uhf
         INTEGER ISPINS,ISPN,ierr,SYMLZ(1000)!,IDI,IDJ,IDK,IDL
         integer OCC(nIrreps), CLOSED(nIrreps), FROZEN(nIrreps)
         INTEGER TMatSize,IUHF
         integer(int64) :: UMatSize
         INTEGER , ALLOCATABLE :: CacheInd(:)
         character(len=*), parameter :: t_r='READFCIINT'
         real(dp) :: diff
         logical :: tbad,tRel
         integer(int64) :: start_ind, end_ind
         integer(int64), parameter :: chunk_size = 1000000
         integer:: bytecount
         NAMELIST /FCI/ NORB,NELEC,MS2,ORBSYM,OCC,CLOSED,FROZEN,&
              ISYM,IUHF,UHF,TREL,SYML,SYMLZ,PROPBITLEN,NPROP

#ifdef _MOLCAS_
         logical :: tExists     !test for existence of input file.
         integer :: isfreeunit  !function returning integer for free unit
#endif

         LWRITE=.FALSE.
         UHF=.FALSE.
         TREL=.false.
         IUHF=0
         PROPBITLEN = 0
         NPROP = 0
         SYMLZ(:) = 0
         ZeroedInt=0
         LzDisallowed=0
         NonZeroInt=0
         iunit = 0

         IF(iProcIndex.eq.0) THEN
#ifdef _MOLCAS_
           call f_Inquire('FCIDMP',tExists)
           if(tExists) then
             iUnit=17
             iUnit=IsFreeUnit(iUnit)
             Call Molcas_Open(iunit,'FCIDMP')
             Rewind(iunit)
             READ(iunit,FCI)
!             write(6,*) 'FCI NAMELIST print 3'
!             WRITE(6,FCI)
             CALL neci_flush(6)
           else
              call Stop_All('InitFromFCID','FCIDUMP file does not exist')
           end if
#else
             iunit = get_free_unit()
             OPEN(iunit,FILE=FCIDUMP_name,STATUS='OLD')
             READ(iunit,FCI)
#endif
         ENDIF
!Now broadcast these values to the other processors (the values are only read in on root)
         CALL MPIBCast(NORB,1)
         CALL MPIBCast(NELEC,1)
         CALL MPIBCast(MS2,1)
         CALL MPIBCast(ORBSYM,1000)
         CALL MPIBCast(SYML,1000)
         CALL MPIBCast(SYMLZ,1000)
         CALL MPIBCast(ISYM,1)
         CALL MPIBCast(IUHF,1)
         CALL MPIBCast(UHF)
         CALL MPIBCast(TREL)
         CALL MPIBCast(PROPBITLEN,1)
         CALL MPIBCast(NPROP,3)
         ! If PropBitLen has been set then assume we're not using an Abelian
         ! symmetry group which has two cycle generators (ie the group has
         ! complex representations).
         TwoCycleSymGens = PropBitLen == 0
         tReltvy = tRel

         ISPINS=2
         IF(UHF.and.(.not.tROHF)) ISPINS=1
         if(tMolpro.and.tUHF) ISPINS=1
         if(tReltvy) ISPINS=1

         IF(tROHF.and.(.not.tMolpro)) THEN
!We are reading in the symmetry of the orbitals in spin-orbitals - we need to change this to spatial orbitals
!NORB is the length in spatial orbitals, LEN is spin orbitals (both spin for UHF)
             NORB=NORB/2
         ENDIF

         IF(tCacheFCIDUMPInts) THEN
!We need to fill the cache. We do this by filling it contiguously, then ordering.
             ALLOCATE(CacheInd(nPairs),stat=ierr)
             CacheInd(:)=1
         ENDIF

         IF(iProcIndex.eq.0) THEN

             write(6,'("Two-electron integrals with a magnitude over ", &
                       &g16.7," are screened")') UMatEps

             if(tMolpro.and.tUHF) then
                 !In molpro, UHF FCIDUMPs are written out as:
                 !1: aaaa
                 !2: bbbb
                 !3: aabb
                 !4: aa
                 !5: bb
                 !with delimiters of 0.000000 0 0 0 0
                 iSpinType=1
             else
                 iSpinType=0
             endif
             TMAT2D(:,:)=(0.0_dp)
             ! Can't use * as need to be backward compatible with existing
             ! FCIDUMP files, some of which have more than 100 basis
             ! functions and so the integer labels run into each other.
             ! This means it won't work with more than 999 basis
             ! functions...

#ifdef __CMPLX
101          READ(iunit,*,END=199) Z,I,J,K,L
#else
101          CONTINUE
             !It is possible that the FCIDUMP can be written out in complex notation, but still only
             !have real orbitals. This occurs with solid systems, where all kpoints are at the
             !gamma point or BZ boundary and have been appropriately rotated. In this case, all imaginary
             !components should be zero to numerical precision.
             if(tRotatedOrbsReal) then
                 !Still need to read in integrals as complex numbers - this FCIDUMP will be from VASP.
                 READ(iunit,*,END=199) CompInt,I,J,K,L
                 Z=real(CompInt,dp)
                 if(abs(aimag(CompInt)).gt.1.0e-7_dp) then
                     call stop_all("READFCIINT","Real orbitals indicated, but imaginary part of integrals larger than 1.0e-7_dp")
                 endif
             else
                 if(tMolpro.or.tReadFreeFormat) then
                     read(iunit,*,END=199) Z,I,J,K,L
                 else
                     READ(iunit,'(1X,G20.12,4I3)',END=199) Z,I,J,K,L
                 endif
             endif
#endif

             ! Remove integrals that are too small
             if (abs(Z) < UMatEps) then
                 if (ZeroedInt < 100) then
                     write(6,'(a,2i4,a,2i4,a)', advance='no') &
                         'Ignoring integral (chem. notation) (', i, j, '|', k, &
                        l, '): '
                     write(6,*) Z
                 else if (ZeroedInt == 100) then
                     write(6,*) 'Ignored more than 100 integrals.'
                     write(6,*) 'Further threshold truncations not reported explicitly'
                 end if
                 ZeroedInt = ZeroedInt + 1
                 goto 101
             end if

             ! If we are fixing Lz symmetry, test if symmetry-zero elements
             ! are being included
             if (tFixLz) then
                 tbad = .false.
                 if (i /= 0 .and. j /= 0 .and. k /= 0 .and. l /= 0) then
                     if (SymLz(i) + symLz(k) /= SymLz(j) + SymLz(l))  &
                         tbad = .true.
                 end if
                 if (i /= 0 .and. j /= 0 .and. k == 0 .and. l == 0) then
                     if (SymLz(i) /= SymLz(j))  &
                         tbad = .true.
                 end if
                 if (tbad) then
                     if (LzDisallowed < 100) then
                         write(6,'(a,2i4,a,2i4,a)', advance='no') &
                             'Ignoring Lz disallowed integral (chem. notation)&
                             & (', i, j, '|', k, l, '): '
                         write(6,*) Z
                     else if (LzDisallowed == 100) then
                         write(6,*) 'Ignored more than 100 integrals.'
                         write(6,*) 'Further threshold truncations not reported explicitly'
                     end if
                     LzDisallowed = LzDisallowed + 1
                     goto 101
                 end if
             end if


             IF(tROHF.and.(.not.tMolpro)) THEN
!The FCIDUMP file is in spin-orbitals - we need to transfer them to spatial orbitals (unless molpro).
                IF(I.ne.0) THEN
                    I = GTID(I)
                ENDIF
                IF(J.ne.0) THEN
                    J = GTID(J)
                ENDIF
                IF(K.ne.0) THEN
                    K = GTID(K)
                ENDIF
                IF(L.ne.0) THEN
                    L = GTID(L)
                ENDIF
             elseif(tMolpro.and.tUHF) then
                if(i.ne.0) then
                    !Need to transfer the orbital index into spin orbital notation
                    if((iSpinType.eq.1).or.(iSpinType.eq.4)) then
                        !aaaa/aa00
                        I=I*2-1
                        J=J*2-1
                        if(iSpinType.eq.1) K=K*2-1   !(just so it doesn't give -1!)
                        if(iSpinType.eq.1) L=L*2-1
                    elseif((iSpinType.eq.2).or.(iSpinType.eq.5)) then
                        !bbbb/bb00
                        I=I*2
                        J=J*2
                        K=K*2
                        L=L*2
                    elseif(iSpinType.eq.3) then
                        !aabb spin type (remember its still in chemical notation!)
                        I=I*2-1
                        J=J*2-1
                        K=K*2
                        L=L*2
                    endif
                endif
             ENDIF
!.. Each orbital in the file corresponds to alpha and beta spinorbitalsa
             IF(I.EQ.0) THEN
                 if(tMolpro.and.tUHF.and.(iSpinType.ne.6)) then
                     if(abs(real(z,dp)).gt.1.0e-8_dp) then
                         call stop_all(t_r,"This is not a delimiter in the FCIDUMP reading")
                     endif
                     iSpinType=iSpinType+1
                 else
!.. Core energy
                     ECORE=real(Z,dp)
                 endif
             ELSEIF(J.EQ.0) THEN
!C.. HF Eigenvalues
!                ARR(I*2-1,2)=real(Z,dp)
!                ARR(I*2,2)=real(Z,dp)
!                ARR(BRR(I*2-1),1)=real(Z,dp)
!                ARR(BRR(I*2),1)=real(Z,dp)
!                LWRITE=.TRUE.
             ELSEIF(K.EQ.0) THEN
!.. 1-e integrals
!.. These are stored as spinorbitals (with elements between different spins being 0
                DO ISPN=1,ISPINS

                    ! Have read in T_ij.  Check it's consistent with T_ji
                    ! (if T_ji has been read in).
                   diff = abs(TMAT2D(ISPINS*I-ISPN+1,ISPINS*J-ISPN+1)-Z)
                   IF(.not. near_zero(TMAT2D(ISPINS*I-ISPN+1,ISPINS*J-ISPN+1)) .and. diff > 1.0e-7_dp) then
                        WRITE(6,*) i,j,Z,TMAT2D(ISPINS*I-ISPN+1,ISPINS*J-ISPN+1)
                        CALL Stop_All("ReadFCIInt","Error filling TMAT - different values for same orbitals")
                   ENDIF

                   TMAT2D(ISPINS*I-ISPN+1,ISPINS*J-ISPN+1)=Z
#ifdef __CMPLX
                   TMAT2D(ISPINS*J-ISPN+1,ISPINS*I-ISPN+1)=conjg(Z)
#else
                   TMAT2D(ISPINS*J-ISPN+1,ISPINS*I-ISPN+1)=Z
#endif
                enddo
             ELSE
!.. 2-e integrals
!.. UMAT is stored as just spatial orbitals (not spinorbitals)
!..  we're reading in (IJ|KL), but we store <..|..> which is <IK|JL>
#ifdef __CMPLX
                Z = UMatConj(I,K,J,L,Z)
#endif
                IF(TUMAT2D) THEN
                    IF(I.eq.J.and.I.eq.K.and.I.eq.L) THEN
                        !<ii|ii>
                        UMAT2D(I,I)=Z
                    ELSEIF((I.eq.J.and.K.eq.L)) THEN
                        !<ij|ij> - coulomb - 1st arg > 2nd arg
                        X=MAX(I,K)
                        Y=MIN(I,K)
                        UMAT2D(Y,X)=Z
                    ELSEIF(I.eq.L.and.J.eq.K) THEN
                        !<ij|ji> - exchange - 1st arg < 2nd arg
                        X=MIN(I,J)
                        Y=MAX(I,J)
                        UMAT2D(Y,X)=Z
                    ELSEIF(I.eq.K.and.J.eq.L) THEN
                        !<ii|jj> - equivalent exchange for real orbs
                        X=MIN(I,J)
                        Y=MAX(I,J)
                        UMAT2D(Y,X)=Z
                    ELSEIF(tRIIntegrals) THEN
                        CALL Stop_All("ReadFCIINTS","we should not be " &
     &                        //"reading in generic 2e integrals from " &
     &                        //"the FCIDUMP file with ri.")
                    ELSEIF(tCacheFCIDUMPInts) THEN
!Read in the FCIDUMP integrals to a cache.
                        IF(tReadFreezeInts) THEN
!Here, we only want to cache the <ij|kj> integrals, since they are the only integrals which are needed for the freezing process.
                          CALL Stop_All("READFCIINTS","Freezing is not yet compatible with caching the FCIDUMP file")
                        ELSE
                         CALL CacheFCIDUMP(I,K,J,L,Z,CacheInd,ZeroedInt,NonZeroInt)

                        ENDIF
                    ELSE
                        NonZeroInt=NonZeroInt+1
!Read in all integrals as normal.
                        UMAT(UMatInd(I,K,J,L))=Z
                    ENDIF
                ELSEIF(tCacheFCIDUMPInts.or.tRIIntegrals) THEN
                    CALL Stop_All("ReadFCIInts","TUMAT2D should be set")
                ELSE
                    UMAT(UMatInd(I,K,J,L))=Z
                    NonZeroInt=NonZeroInt+1
                ENDIF
             ENDIF
!             IF(I.NE.0) GOTO 101
             GOTO 101
199          CONTINUE

             CLOSE(iunit)
             if(tMolpro.and.tUHF.and.(iSpinType.ne.6)) then
                 call stop_all(t_r,"Error in reading UHF FCIDUMP from molpro")
             endif
         ENDIF

!Now broadcast the data read in
         CALL MPIBCast(ZeroedInt,1)
         CALL MPIBCast(LzDisallowed,1)
         CALL MPIBCast(NonZeroInt,1)
         CALL MPIBCast(ECore,1)
!Need to find out size of TMAT before we can BCast
         CALL CalcTMATSize(nBasis,TMATSize)
         CALL MPIBCast(TMAT2D,TMATSize)
         IF(TUMAT2D) THEN
!Broadcast TUMAT2D...
             CALL MPIBCast(UMAT2D,nStates**2)
         ENDIF
         IF((.not.tRIIntegrals).and.(.not.tCacheFCIDUMPInts)) THEN
             CALL GetUMATSize(nBasis,NEl,UMatSize)

             ! If we are on a 64bit system, the maximum dimensions for MPI are
             ! still limited by 32bit limits.
             ! --> We need to loop around this
             start_ind = 1
             end_ind = min(UMatSize, chunk_size)

             do while(start_ind <= UMatSize)
                !use MPI_BYTE for transfer to be independent of the data type of UMat
                bytecount = int(end_ind - start_ind + 1_int64) * int(sizeof(UMat(1)))
                call MPIBCast_inter_byte(UMat(start_ind),bytecount)
                start_ind = end_ind + 1
                end_ind = min(UMatSize, end_ind + chunk_size)
             end do

             !make sure the shared memory data is synchronized on all tasks
             call shared_sync_mpi(umat_win)
         ENDIF

         IF(tCacheFCIDUMPInts) THEN
!Need to broadcast the cache...
             CALL MPIBCast(UMATLABELS,nSlots*nPairs)
             CALL MPIBCast(UMatCacheData,nTypes*nSlots*nPairs)
         ENDIF

         if (ZeroedInt /= 0 .and. iProcIndex == 0) then
             write(6,*) 'Number of removed two-index integrals: ', zeroedint
         end if
         if (LzDisallowed /= 0 .and. iProcIndex == 0) then
             write(6,*) 'Number of Lz disallowed two-index integrals: ', &
                 LzDisallowed
         end if
         write(6,*) 'Number of non-zero integrals: ', NonZeroInt


         IF(tCacheFCIDUMPInts) THEN
             WRITE(6,*) "Ordering cache..."
             CALL FillUpCache()
             DEALLOCATE(CacheInd)
!             CALL DumpUMatCache()
!             do i=1,norb
!                 WRITE(iunit,*) UMAT2D(:,i)
!             enddo
         ENDIF
!.. If we've changed the eigenvalues, we write out the basis again
!         IF(LWRITE) THEN
!            WRITE(6,*) "1-electron energies have been read in."
!            CALL WRITEBASIS(6,G1,NBASIS,ARR,BRR)
!         ENDIF
         RETURN
      END SUBROUTINE READFCIINT


      !This is a copy of the routine above, but now for reading in binary files of integrals
      SUBROUTINE READFCIINTBIN(UMAT,ECORE)
         use constants, only: dp,int64,sizeof_int
         use SystemData, only: Symmetry, BasisFN
         USE UMatCache , only : UMatInd,UMAT2D,TUMAT2D
         use OneEInts, only: TMatind,TMat2D,TMATSYM
         use util_mod, only: get_free_unit
         IMPLICIT NONE
         real(dp), intent(out) :: ECORE
         HElement_t(dp), intent(out) :: UMAT(*)
         HElement_t(dp) Z
         integer(int64) MASK,IND
         INTEGER I,J,K,L,X,Y, iunit
         LOGICAL LWRITE
         character(*), parameter :: this_routine = 'READFCIINTBIN'
         !NAMELIST /FCI/ NORB,NELEC,MS2,ORBSYM,ISYM,UHF
         LWRITE=.FALSE.
         !UHF=.FALSE.
         iunit = get_free_unit()
         !OPEN(iunit,FILE='FCISYM',STATUS='OLD',FORM='FORMATTED')
         !READ(iunit,FCI)
         !CLOSE(iunit)
         OPEN(iunit,FILE=FCIDUMP_name,STATUS='OLD',FORM='UNFORMATTED')


         MASK=(2**16)-1
         !IND contains all the indices in an integer(int64) - use mask of 16bit to extract them
101      READ(iunit,END=199) Z,IND
         L=int(iand(IND,MASK),sizeof_int)
         IND=Ishft(IND,-16)
         K=int(iand(IND,MASK),sizeof_int)
         IND=Ishft(IND,-16)
         J=int(iand(IND,MASK),sizeof_int)
         IND=Ishft(IND,-16)
         I=int(iand(IND,MASK),sizeof_int)
!         I=Index(I)
!         J=Index(J)
!         K=Index(K)
!         L=Index(L)

!.. Each orbital in the file corresponds to alpha and beta spinorbitalsa
         IF(I.EQ.0) THEN
!.. Core energy
            ECORE=real(Z,dp)
         ELSEIF(J.EQ.0) THEN
!C.. HF Eigenvalues
!            ARR(I*2-1,2)=real(Z,dp)
!            ARR(I*2,2)=real(Z,dp)
!            ARR(BRR(I*2-1),1)=real(Z,dp)
!            ARR(BRR(I*2),1)=real(Z,dp)
!            LWRITE=.TRUE.
         ELSEIF(K.EQ.0) THEN
!.. 1-e integrals
!.. These are stored as spinorbitals (with elements between different spins being 0
            TMAT2D(2*I-1,2*J-1)=Z
            TMAT2D(2*I,2*J)=Z
            TMAT2D(2*J-1,2*I-1)=Z
            TMAT2D(2*J,2*I)=Z
         ELSE
!.. 2-e integrals
!.. UMAT is stored as just spatial orbitals (not spinorbitals)
!..  we're reading in (IJ|KL), but we store <..|..> which is <IK|JL>
            UMAT(UMatInd(I,K,J,L))=Z
         ENDIF
!         WRITE(14,'(1X,F20.12,4I3)') Z,I,J,K,L
         IF(I.NE.0) GOTO 101
199      CONTINUE
         CLOSE(iunit)
! If we've changed the eigenvalues, we write out the basis again
!         IF(LWRITE) THEN
!            WRITE(6,*) "1-electron energies have been read in."
!            CALL WRITEBASIS(6,G1,NBASIS,ARR,BRR)
!         ENDIF
         RETURN
      END SUBROUTINE READFCIINTBIN

      SUBROUTINE ReadPropInts(iProp,nBasis,iNumProp,PropFile,CoreVal,OneElInts)

      use constants, only: dp, int64, iout
      use util_mod, only: get_free_unit
      use SymData, only: PropBitLen,nProp
      use SystemData, only: UMatEps, tROHF, tReltvy
      use Parallel_neci, only : iProcIndex,MPIBcast

      implicit none
      integer, intent(in) :: iProp, nBasis
      HElement_t(dp) :: OneElInts(nBasis,nBasis)
      HElement_t(dp) z
      real(dp) :: CoreVal
      integer :: i,j,k,l,iunit, iNumProp
      integer :: NORB,NELEC,MS2,ISYM,SYML(1000),IUHF
      integer(int64) :: ORBSYM(1000)
      integer :: iSpins,ispn,SYMLZ(1000)
      integer(int64) :: ZeroedInt
      integer :: IntSize
      real(dp) :: diff, core
      character(len=100) :: file_name, PropFile
      logical :: TREL,UHF
      character(*), parameter :: t_r='ReadPropInts'
      NAMELIST /FCI/ NORB,NELEC,MS2,ORBSYM,ISYM,IUHF,UHF,TREL,SYML,SYMLZ,PROPBITLEN,NPROP

      ZeroedInt = 0
      UHF = .false.

      if(iProcIndex.eq.0) then
          iunit = get_free_unit()
          file_name = PropFile
          write(iout,*) 'Reading integral from the file:', trim(file_name)
          open(iunit,FILE=file_name,STATUS='OLD')
          read(iunit,FCI)
      end if

!Now broadcast these values to the other processors (the values are only read in on root)
      CALL MPIBCast(NORB,1)
      CALL MPIBCast(NELEC,1)
      CALL MPIBCast(MS2,1)
      CALL MPIBCast(ORBSYM,1000)
      CALL MPIBCast(SYML,1000)
      CALL MPIBCast(SYMLZ,1000)
      CALL MPIBCast(ISYM,1)
      CALL MPIBCast(IUHF,1)
      CALL MPIBCast(UHF)
      CALL MPIBCast(TREL)
      CALL MPIBCast(PROPBITLEN,1)
      CALL MPIBCast(NPROP,3)

      core = 0.0d0
      iSpins=2
      IF((UHF.and.(.not.tROHF)).or.tReltvy) ISPINS=1

      if(iProcIndex.eq.0) then
101       continue
          read(iunit,*,END=199) z,i,j,k,l

          ! Remove integrals that are too small
          if (abs(z) < UMatEps) then
              if (ZeroedInt < 100) then
                  write(6,'(a,2i4,a,2i4,a)', advance='no') &
                      'Ignoring integral (chem. notation) (', i, j, '|', k, &
                     l, '): '
                  write(6,*) z
              else if (ZeroedInt == 100) then
                  write(6,*) 'Ignored more than 100 integrals.'
                  write(6,*) 'Further threshold truncations not reported explicitly'
              end if
              ZeroedInt = ZeroedInt + 1
              goto 101
          end if

          if (i.eq.0) then
! Reading the zero-electron part of the integrals
              core=real(z,dp)
          elseif (k.eq.0) then
! Reading the one-electron part of the integrals
              do ispn=1,iSpins
                  ! Have read in T_ij.  Check it's consistent with T_ji
                  ! (if T_ji has been read in).
                  diff = abs(OneElInts(iSpins*i-ispn+1,iSpins*j-ispn+1)-z)
                  if(abs(OneElInts(iSpins*i-ispn+1,iSpins*j-ispn+1)).gt. 0.0d-6 .and. diff > 1.0e-7_dp) then
                  write(6,*) i,j,z,OneElInts(iSpins*i-ispn+1,iSpins*j-ispn+1)
                  call Stop_All(t_R,"Error filling OneElInts - different values for same orbitals")
                  endif

                  OneElInts(iSpins*I-ispn+1,iSpins*J-ispn+1)=z
#ifdef __CMPLX
                  OneElInts(iSpins*J-ispn+1,iSpins*I-ispn+1)=conjg(z)
#else
                  OneElInts(iSpins*J-ispn+1,iSpins*I-ispn+1)=z
#endif
              enddo

          else
              call stop_all(t_r,'Cannot currently deal with two-electron properties')
          endif
          goto 101
      endif
199       continue
      if(iProcIndex.eq.0) close(iunit)

      CoreVal = core

      END SUBROUTINE ReadPropInts

end module read_fci
