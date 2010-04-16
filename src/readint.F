
      SUBROUTINE INITFROMFCID(NEL,NBASISMAX,LEN,LMS,TBIN)
         use SystemData , only : tNoSymGenRandExcits,lNoSymmetry,tROHF
         use SystemData , only : tStoreSpinOrbs
         use Parallel
         IMPLICIT NONE
         INTEGER NEL,nBasisMax(5,*),LEN,LMS,SYMLZ(1000)
         INTEGER NORB,NELEC,MS2,ORBSYM(1000),ISYM,i,SYML(1000)
         LOGICAL TBIN,exists,UHF
         CHARACTER*3 :: fmat
         NAMELIST /FCI/ NORB,NELEC,MS2,ORBSYM,ISYM,UHF,SYML,SYMLZ
         UHF=.FALSE.
         fmat='NO'
         IF(iProcIndex.eq.0) THEN
             IF(TBIN) THEN
                INQUIRE(FILE='FCISYM',EXIST=exists)
                IF(.not.exists) THEN
                    CALL Stop_All('InitFromFCID',                       &
     &                   'FCISYM file does not exist')
                ENDIF
                INQUIRE(FILE='FCIDUMP',EXIST=exists,FORMATTED=fmat)
                IF(.not.exists) THEN
                    CALL Stop_All('INITFROMFCID',                       &
     &                    'FCIDUMP file does not exist')
                ENDIF
!                IF(fmat=='YES') THEN
!                    STOP 'FCIDUMP is not unformatted, but TBIN true'
!                ENDIF
                OPEN(8,FILE='FCISYM',STATUS='OLD',FORM='FORMATTED')
                READ(8,FCI)
             ELSE
                INQUIRE(FILE='FCIDUMP',EXIST=exists,UNFORMATTED=fmat)
                IF(.not.exists) THEN
                    CALL Stop_All('InitFromFCID',                       &
     &                  'FCIDUMP file does not exist')
                ENDIF
!                IF(fmat=='YES') THEN
!                    STOP 'FCIDUMP is not formatted, but TBIN false'
!                ENDIF
                OPEN(8,FILE='FCIDUMP',STATUS='OLD',FORM='FORMATTED')
                READ(8,FCI)
             ENDIF
             CLOSE(8)
         ENDIF

!Now broadcast these values to the other processors
         CALL MPIIBCast_Scal(NORB,0)
         CALL MPIIBCast_Scal(NELEC,0)
         CALL MPIIBCast_Scal(MS2,0)
         CALL MPIIBCast(ORBSYM,1000,0)
         CALL MPIIBCast(SYML,1000,0)
         CALL MPIIBCast(SYMLZ,1000,0)
         CALL MPIIBCast_Scal(ISYM,0)
         CALL MPILBCast_Scal(UHF,0)

         IF(tROHF.and.(.not.UHF)) THEN
             CALL Stop_All("INITFROMFCID","ROHF specified, but FCIDUMP "&
     &         //"is not in a high-spin format.")
         ENDIF


         DO i=1,NORB
             IF(ORBSYM(i).eq.0) THEN
                 WRITE(6,*) "** WARNING **"
                 WRITE(6,*) "** Unconverged symmetry of orbitals **"
                 WRITE(6,*) "** Turning symmetry off for rest of run **"
                 tNoSymGenRandExcits=.true.
                 lNoSymmetry=.true.
                 EXIT
             ENDIF
         ENDDO
         IF(NELEC.NE.NEL) THEN
             WRITE(6,*)                                                 &
     &      '*** WARNING: NEL in FCIDUMP differs from input file ***'   
            WRITE(6,*) ' NUMBER OF ELECTRONS : ' , NEL
         ENDIF
C         NEL=NELEC
         IF(LMS.NE.MS2) THEN   
            WRITE(6,*)                                                  &
     &      '*** WARNING: LMS in FCIDUMP differs from input file ***'   
            LMS=MS2                                                     
            WRITE(6,*) ' BASIS MS : ' , LMS
         ENDIF
         IF(UHF) then
            LEN=NORB
         ELSE
            LEN=2*NORB
         ENDIF
         NBASISMAX(1:5,1:3)=0
C         NBASISMAX(1,1)=1
C         NBASISMAX(1,2)=NORB
         NBASISMAX(1,1)=0
         NBASISMAX(1,2)=0
         NBASISMAX(1,3)=2
         NBASISMAX(4,1)=-1
         NBASISMAX(4,2)=1
         
!the indicator of UHF (which becomes ISpinSkip or ISS later
         IF(UHF.and.(.not.tROHF)) then
            NBASISMAX(2,3)=1
            tStoreSpinOrbs=.true.   !indicate that we are storing the orbitals in umat as spin-orbitals
         else
            NBASISMAX(2,3)=2
            tStoreSpinOrbs=.false.   !indicate that we are storing the orbitals in umat as spatial-orbitals
         endif

C.. SHow that there's no momentum conservation
         NBASISMAX(3,3)=1
         RETURN
      END


      SUBROUTINE GETFCIBASIS(NBASISMAX,ARR,BRR,G1,LEN,TBIN)
         use SystemData, only: BasisFN,BasisFNSize,Symmetry,NullBasisFn
         use SystemData, only: tCacheFCIDUMPInts,tROHF,tFixLz,iMaxLz
         use UMatCache, only: nSlotsInit,CalcNSlotsInit
         use UMatCache, only: GetCacheIndexStates,GTID
         use Parallel
         IMPLICIT NONE
         INTEGER nBasisMax(5,*),BRR(LEN),LEN
         TYPE(BasisFN) G1(LEN)
         REAL*8 ARR(LEN,2)
         REAL*8 Z
         INTEGER*8 IND,MASK
         INTEGER I,J,K,L,I1
         INTEGER ISYMNUM,ISNMAX,SYMMAX,SYMLZ(1000)
         INTEGER NORB,NELEC,MS2,ORBSYM(1000),ISYM,ISPINS,ISPN,SYML(1000)
         INTEGER Counter(1:8),nPairs,iErr,MaxnSlot,MaxIndex
         INTEGER , ALLOCATABLE :: MaxSlots(:)
         LOGICAL TBIN,UHF
         NAMELIST /FCI/ NORB,NELEC,MS2,ORBSYM,ISYM,UHF,SYML,SYMLZ
         UHF=.FALSE.
         IF(iProcIndex.eq.0) THEN
             IF(TBIN) THEN
                OPEN(8,FILE='FCISYM',STATUS='OLD',FORM='FORMATTED')
                READ(8,FCI)
                CLOSE(8)
                OPEN(8,FILE='FCIDUMP',STATUS='OLD',FORM='UNFORMATTED')
             ELSE
                OPEN(8,FILE='FCIDUMP',STATUS='OLD',FORM='FORMATTED')
                READ(8,FCI)
             ENDIF
         ENDIF

!Now broadcast these values to the other processors (the values are only read in on root)
         CALL MPIIBCast_Scal(NORB,0)
         CALL MPIIBCast_Scal(NELEC,0)
         CALL MPIIBCast_Scal(MS2,0)
         CALL MPIIBCast(ORBSYM,1000,0)
         CALL MPIIBCast(SYML,1000,0)
         CALL MPIIBCast(SYMLZ,1000,0)
         CALL MPIIBCast_Scal(ISYM,0)
         CALL MPILBCast_Scal(UHF,0)

         ISYMNUM=0
         ISNMAX=0
         ISPINS=2
         IF(UHF.and.(.not.tROHF)) ISPINS=1

         IF(tROHF) THEN
             WRITE(6,*) "Reading in Spin-orbital FCIDUMP but storing as"&
     &          //" spatial orbitals."
!             WRITE(6,*) "*** Warning - Fock energy of orbitals will be "&
!     &                //"incorrect ***"
!We are reading in the symmetry of the orbitals in spin-orbitals - we need to change this to spatial orbitals
             DO i=1,NORB-1,2
                 IF(ORBSYM(i).ne.ORBSYM(i+1)) THEN
                     CALL Stop_All("GetFCIBASIS","Spin-orbitals are not"&
     &                  //" ordered in symmetry pairs")
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
         ENDIF
        
         IF(LEN.NE.ISPINS*NORB) STOP 'LEN .NE. NORB in GETFCIBASIS'
         G1(1:LEN)=NullBasisFn
         ARR=0.d0

!!Calculate a key, which will reorder the orbitals in symmetry order.
!!Counter is the cumulative number of symmetries, just to indicate the allowed current index of an orbital
!         Counter(:)=0
!         SymMax=1
!         do i=1,norb
!             IF(ORBSYM(i).gt.SymMax) SymMax=ORBSYM(i)
!             Counter(ORBSYM(i))=Counter(ORBSYM(i))+1
!         enddo
!!         WRITE(6,*) Counter(:)
!         do i=SymMax,2,-1
!             Counter(i)=1
!             do j=1,i-1
!                 Counter(i)=Counter(i)+Counter(j)
!             enddo
!         enddo
!         Counter(1)=1
!!Check that all orbitals are accounted for.
!         j=0
!         do i=1,norb
!             IF(ORBSYM(i).eq.SymMax) j=j+1
!         enddo
!         IF((Counter(SymMax)+j-1).ne.norb) THEN
!             WRITE(6,*) "***",SymMax,norb,j
!             WRITE(6,*) Counter(1:SymMax)
!             CALL Stop_All("GETFCIBASIS","Counter incorrectly set up")
!         ENDIF
!!Now work out what the new indices are to be called
!         Index(1:norb)=0
!         do i=1,norb
!             Index(i)=Counter(OrbSym(i))
!             Counter(OrbSym(i))=Counter(OrbSym(i))+1
!         enddo
!!Now need to sort OrbSym so that the symmetries correspond to the correct orbitals.
!         CALL NECI_SORTI(norb,OrbSym(1:norb))

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
             CALL FLUSH(6)
             ALLOCATE(MaxSlots(nPairs),stat=ierr)
             MaxSlots(:)=0
         ENDIF


         IF(iProcIndex.eq.0) THEN
!Just read in integrals on head node.

             IF(TBIN) THEN
                IF(UHF.or.tROHF) STOP 'UHF Bin read not functional'
                MASK=(2**16)-1
                
                !IND contains all the indices in an integer*8 - use mask of 16bit to extract them
2               READ(8,END=99) Z,IND
                L=iand(IND,MASK)
                IND=Ishft(IND,-16)
                K=iand(IND,MASK)
                IND=Ishft(IND,-16)
                J=iand(IND,MASK)
                IND=Ishft(IND,-16)
                I=iand(IND,MASK)

!                I=Index(I)
!                J=Index(J)
!                K=Index(K)
!                L=Index(L)
            
C.. Each orbital in the file corresponds to alpha and beta spinorbitals
             !Fill ARR with the energy levels
                IF(I.NE.0.AND.K.EQ.0.AND.I.EQ.J) THEN
                    IF(I.GT.1) THEN
                        IF(ORBSYM(I).NE.ORBSYM(I-1)) THEN
                            IF(ISYMNUM.GT.ISNMAX) ISNMAX=ISYMNUM
                            ISYMNUM=0
                        ENDIF
                    ENDIF
                    ISYMNUM=ISYMNUM+1
                    ARR(2*I-1,1)=Z
                    ARR(2*I,1)=Z
                ELSEIF(I.NE.0.AND.K.EQ.0.AND.J.EQ.0) THEN
                    ARR(2*I-1,1)=Z
                    ARR(2*I,1)=Z
                ENDIF
C.. At the moment we're ignoring the core energy
                IF(I.NE.0) GOTO 2

             ELSE   !Reading in formatted FCIDUMP file

1               READ(8,'(1X,G20.12,4I3)',END=99) Z,I,J,K,L

                IF(tROHF) THEN
!The FCIDUMP file is in spin-orbitals - we need to transfer them to spatial orbitals.
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

!                I=Index(I)
!                J=Index(J)
!                K=Index(K)
!                L=Index(L)
C.. Each orbital in the file corresponds to alpha and beta spinorbitals
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
                    IF(ISPINS.eq.1) THEN
                        ARR(I1,1)=Z
                    ELSEIF(tROHF) THEN
                        ARR(I,1)=Z
                    ELSE
                        ARR(2*I1,1)=Z
                        ARR(2*I1-1,1)=Z
                    ENDIF
!                    DO ISPN=1,ISPINS
!                        ARR(ISPINS*I1-ISPN+1,1)=Z
!                    ENDDO

                ELSEIF(I1.NE.0.AND.K.EQ.0.AND.J.EQ.0) THEN
!                    WRITE(6,*) I
                    IF(ISPINS.eq.1) THEN
                        ARR(I1,1)=Z
                    ELSEIF(tROHF) THEN
                        ARR(I,1)=Z
                    ELSE
                        ARR(2*I1,1)=Z
                        ARR(2*I1-1,1)=Z
                    ENDIF

!                    DO ISPN=1,ISPINS
!                        ARR(ISPINS*I-ISPN+1,1)=Z
!                    ENDDO

                ENDIF
C.. At the moment we're ignoring the core energy
                IF(I1.NE.0) GOTO 1
             ENDIF
99           CONTINUE
             CLOSE(8)

         ENDIF

!We now need to broadcast all the information we've just read in...
         IF(tCacheFCIDUMPInts) THEN
             CALL MPIIBCast(MaxSlots,nPairs,0)
         ENDIF
         CALL MPIIBCast_Scal(ISNMAX,0)
         CALL MPIIBCast_Scal(ISYMNUM,0)
         CALL MPIDBCastArr(Arr,LEN*2,0)
         
         SYMMAX=1
         iMaxLz=0
         DO I=1,NORB
            DO ISPN=1,ISPINS
                BRR(ISPINS*I-ISPN+1)=ISPINS*I-ISPN+1
                G1(ISPINS*I-ISPN+1)%Sym%s=ORBSYM(I)-1
                IF(tFixLz) THEN
                    G1(ISPINS*I-ISPN+1)%Ml=SYMLZ(I)
                ELSE
                    G1(ISPINS*I-ISPN+1)%Ml=0
                ENDIF
C.. set momentum to 0
                G1(ISPINS*I-ISPN+1)%k(1)=0
                G1(ISPINS*I-ISPN+1)%k(2)=0
                G1(ISPINS*I-ISPN+1)%k(3)=0
                G1(ISPINS*I-ISPN+1)%Ms=-MOD(ISPINS*I-ISPN+1,2)*2+1
                IF(SYMMAX.lt.ORBSYM(I)) SYMMAX=ORBSYM(I)
                IF(abs(SYMLZ(I)).gt.iMaxLz) iMaxLz=abs(SYMLZ(I))
            ENDDO
         ENDDO
         IF(.not.tFixLz) iMaxLz=0
         IF(SYMMAX.EQ.1) THEN
         ELSEIF(SYMMAX.LE.2) THEN
         ELSEIF(SYMMAX.LE.4) THEN
            SYMMAX=4
         ELSEIF(SYMMAX.LE.8) THEN
            SYMMAX=8
         ENDIF
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
             WRITE(6,*) "Maximum number of slots needed in cache is: ", &
     &           nSlotsInit
             WRITE(6,*) "Index corresponding to this maximum number:",  &
     &           MaxIndex
             CALL GetCacheIndexStates(MaxIndex,i,k)
             WRITE(6,"(A,2I6)") "This corresponds to an (i,k) pair of ",&
     &           i,k
         ENDIF
!         WRITE(6,*) Arr(:,1)
         RETURN
      END

!tReadFreezeInts only matters when we are cacheing the FCIDUMP file.
!It is set if we want to cache the integrals to enable the freezing routine to take place, i.e. the <ij|kj> integrals.
!The UMAT2D integrals will also be read in in this case.
!If tReadFreezeInts is false, then if we are cacheing the FCIDUMP file, then we will read and cache all the integrals.
      SUBROUTINE READFCIINT(UMAT,NBASIS,ECORE,ARR,BRR,G1,               &
     &   tReadFreezeInts)
         use constants, only: dp
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB,NEl
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SystemData, only: UMatEps,tUMatEps,tCacheFCIDUMPInts
         use SystemData, only: tRIIntegrals,nBasisMax,tROHF
         USE UMatCache, only: UMatInd,UMAT2D,TUMAT2D,nPairs,CacheFCIDUMP
         USE UMatCache, only: FillUpCache,GTID,nStates,nSlots,nTypes
         USE UMatCache, only: UMatCacheData,UMatLabels,GetUMatSize
         use OneEInts, only: TMatind,TMat2D,TMATSYM,TSTARSTORE
         use OneEInts, only: CalcTMatSize
         use Parallel
         IMPLICIT NONE
         INTEGER NBASIS,ZeroedInt,NonZeroInt
         REAL*8 ECORE,ARR(NBASIS,2)
         REAL*8 UMAT(*)
         REAL*8 Z
         HElement_t UMatEl
         TYPE(BasisFN) G1(*)
         INTEGER I,J,K,L,BRR(NBASIS),X,Y,A,B,iCache,iCacheI,iType
         INTEGER NORB,NELEC,MS2,ORBSYM(1000),ISYM,SYMMAX,SYML(1000)
         LOGICAL LWRITE,UHF,tAddtoCache,GetCachedUMatEl,tReadFreezeInts
         INTEGER ISPINS,ISPN,ISPN2,ierr,SYMLZ(1000)!,IDI,IDJ,IDK,IDL
         INTEGER Counter(1:8),UMatSize,TMatSize
         INTEGER , ALLOCATABLE :: CacheInd(:)
         NAMELIST /FCI/ NORB,NELEC,MS2,ORBSYM,ISYM,UHF,SYML,SYMLZ
         LWRITE=.FALSE.
         UHF=.FALSE.
         ZeroedInt=0
         NonZeroInt=0
         
         IF(iProcIndex.eq.0) THEN
             OPEN(8,FILE='FCIDUMP',STATUS='OLD')
             READ(8,FCI)
         ENDIF
!Now broadcast these values to the other processors (the values are only read in on root)
         CALL MPIIBCast_Scal(NORB,0)
         CALL MPIIBCast_Scal(NELEC,0)
         CALL MPIIBCast_Scal(MS2,0)
         CALL MPIIBCast(ORBSYM,1000,0)
         CALL MPIIBCast(SYML,1000,0)
         CALL MPIIBCast(SYMLZ,1000,0)
         CALL MPIIBCast_Scal(ISYM,0)
         CALL MPILBCast_Scal(UHF,0)

         ISPINS=2
         IF(UHF.and.(.not.tROHF)) ISPINS=1

         IF(tROHF) THEN
!We are reading in the symmetry of the orbitals in spin-orbitals - we need to change this to spatial orbitals
!NORB is the length in spatial orbitals, LEN is spin orbitals (both spin for UHF)
             NORB=NORB/2
         ENDIF

!!Calculate a key, which will reorder the orbitals in symmetry order.
!!Counter is the cumulative number of symmetries, just to indicate the allowed current index of an orbital
!         Counter(:)=0
!         SymMax=1
!         do i=1,NORB
!             IF(ORBSYM(i).gt.SymMax) SymMax=ORBSYM(i)
!             Counter(ORBSYM(i))=Counter(ORBSYM(i))+1
!         enddo
!         do i=SymMax,2,-1
!             Counter(i)=1
!             do j=1,i-1
!                 Counter(i)=Counter(i)+Counter(j)
!             enddo
!         enddo
!         Counter(1)=1
!!Check that all orbitals are accounted for.
!         j=0
!         do i=1,NORB
!             IF(ORBSYM(i).eq.SymMax) j=j+1
!         enddo
!         IF((Counter(SymMax)+j-1).ne.NORB) THEN
!             CALL Stop_All("GETFCIBASIS","Counter incorrectly set up")
!         ENDIF
!!Now work out what the new indices are to be called
!         Index(1:NORB)=0
!         do i=1,NORB
!             Index(i)=Counter(OrbSym(i))
!             Counter(OrbSym(i))=Counter(OrbSym(i))+1
!         enddo
!!Now need to sort OrbSym so that the symmetries correspond to the correct orbitals.
!         CALL NECI_SORTI(NORB,OrbSym(1:NORB))
                
         IF(tCacheFCIDUMPInts) THEN
!We need to fill the cache. We do this by filling it contiguously, then ordering.
             ALLOCATE(CacheInd(nPairs),stat=ierr)
             CacheInd(:)=1
         ENDIF


         IF(iProcIndex.eq.0) THEN
             IF(.not.TSTARSTORE) THEN
                 TMAT2D(:,:)=(0.D0)
             ENDIF
101          READ(8,'(1X,G20.12,4I3)',END=199) Z,I,J,K,L
             IF(tROHF) THEN
!The FCIDUMP file is in spin-orbitals - we need to transfer them to spatial orbitals.
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
             ENDIF
!             I=Index(I)
!             J=Index(J)
!             K=Index(K)
!             L=Index(L)
C.. Each orbital in the file corresponds to alpha and beta spinorbitalsa
             IF(I.EQ.0) THEN
C.. Core energy
                ECORE=Z
             ELSEIF(J.EQ.0) THEN
!C.. HF Eigenvalues
!                ARR(I*2-1,2)=Z
!                ARR(I*2,2)=Z
!                ARR(BRR(I*2-1),1)=Z
!                ARR(BRR(I*2),1)=Z
!                LWRITE=.TRUE.
             ELSEIF(K.EQ.0) THEN
C.. 1-e integrals
                IF(TSTARSTORE) THEN
! If TSTARSTORE, the one-el integrals are stored in symmetry classes, as spatial orbitals
                    TMATSYM(TMatInd(2*J,2*I))=Z
                ELSE
C.. These are stored as spinorbitals (with elements between different spins being 0
                    DO ISPN=1,ISPINS

                       IF((TMAT2D(ISPINS*I-ISPN+1,ISPINS*J-ISPN+1)      &
     &                      .ne.0.D0).and.abs(TMAT2D(ISPINS*I-ISPN+1,   &
     &                      ISPINS*J-ISPN+1)-Z).gt.1.D-7) THEN
                WRITE(6,*) i,j,Z,TMAT2D(ISPINS*I-ISPN+1,ISPINS*J-ISPN+1)
                           CALL Stop_All("ReadFCIInt","Error filling "  &
     &                    //"TMAT - different values for same orbitals")
                       ENDIF

                       TMAT2D(ISPINS*I-ISPN+1,ISPINS*J-ISPN+1)=Z
                       TMAT2D(ISPINS*J-ISPN+1,ISPINS*I-ISPN+1)=Z
                    enddo
                ENDIF
             ELSE
C.. 2-e integrals
C.. UMAT is stored as just spatial orbitals (not spinorbitals)
C..  we're reading in (IJ|KL), but we store <..|..> which is <IK|JL>
C.. AJWT removed the restriction to TSTARSTORE
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
                          CALL Stop_All("READFCIINTS","Freezing is not "&
     &                 //"yet compatible with caching the FCIDUMP file")
                        ELSE
                         CALL CacheFCIDUMP(I,K,J,L,Z,CacheInd,ZeroedInt,&
     &                          NonZeroInt)

                        ENDIF
                    ELSE
!Read in all integrals as normal.
                        UMAT(UMatInd(I,K,J,L,0,0))=Z
                    ENDIF
                ELSEIF(TSTARSTORE.and.(.not.TUMAT2D)) THEN
                    STOP 'Need UMAT2D with TSTARSTORE'
                ELSEIF(tCacheFCIDUMPInts.or.tRIIntegrals) THEN
                    CALL Stop_All("ReadFCIInts","TUMAT2D should be set")
                ELSE
                    IF(tUMatEps) THEN
!We have an epsilon cutoff for the size of the two-electron integrals - UMatEps
                        IF(abs(Z).lt.UMatEps) THEN
                            UMAT(UMatInd(I,K,J,L,0,0))=0.D0
                            ZeroedInt=ZeroedInt+1
                        ELSE
                            UMAT(UMatInd(I,K,J,L,0,0))=Z
                            NonZeroInt=NonZeroInt+1
                        ENDIF
                    ELSE
                        UMAT(UMatInd(I,K,J,L,0,0))=Z
                    ENDIF
                ENDIF
             ENDIF
             IF(I.NE.0) GOTO 101
199          CONTINUE
             CLOSE(8)
         ENDIF

!Now broadcast the data read in
         CALL MPIIBCast_Scal(ZeroedInt,0)
         CALL MPIIBCast_Scal(NonZeroInt,0)
         CALL MPIDBCast_Scal(ECore,0)
!Need to find out size of TMAT before we can BCast
         CALL CalcTMATSize(nBasis,TMATSize)
         IF(tStarStore) THEN
             CALL MPIHElemBCast(TMATSYM,TMATSize,0)
         ELSE
             CALL MPIHElemBCast(TMAT2D,TMATSize,0)
         ENDIF
         IF(TUMAT2D) THEN
!Broadcast TUMAT2D...
             CALL MPIHElemBCast(UMAT2D,nStates**2,0)
         ENDIF
         IF((.not.tRIIntegrals).and.(.not.tCacheFCIDUMPInts)) THEN
             CALL GetUMATSize(nBasis,NEl,iSpins,UMatSize)
             CALL MPIDBCastArr(UMAT,UMatSize,0)    !This is not an , as it is actually passed in as a real*8, even though it is HElem in IntegralsData
         ENDIF
         IF(tCacheFCIDUMPInts) THEN
!Need to broadcast the cache...
             CALL MPIIBCast(UMATLABELS,nSlots*nPairs,0)
             CALL MPIHElemBCast(UMatCacheData,nTypes*nSlots*nPairs,0)
         ENDIF
             

         IF(tUMatEps) THEN
!Write out statistics from zeroed integrals.
             WRITE(6,"(A,G20.10,A)") "*** Zeroing all two-electron "    &
     &          //"integrals with a magnitude of over ",UMatEps," ***" 
             WRITE(6,*) ZeroedInt+NonZeroInt," 2E integrals read in..."
             WRITE(6,*) ZeroedInt," integrals zeroed..."
             WRITE(6,*) REAL(100*ZeroedInt)/REAL(NonZeroInt+ZeroedInt), &
     %          " percent of 2E integrals zeroed."
         ENDIF
         IF(tCacheFCIDUMPInts) THEN
             WRITE(6,*) "Ordering cache..."
             CALL FillUpCache()
             DEALLOCATE(CacheInd)
!             CALL DumpUMatCache(nBasis,G1)
!             do i=1,norb
!                 WRITE(8,*) UMAT2D(:,i)
!             enddo
         ENDIF
C.. If we've changed the eigenvalues, we write out the basis again
!         IF(LWRITE) THEN
!            WRITE(6,*) "1-electron energies have been read in."
!            CALL WRITEBASIS(6,G1,NBASIS,ARR,BRR)
!         ENDIF
         RETURN
      END


      !This is a copy of the routine above, but now for reading in binary files of integrals
      SUBROUTINE READFCIINTBIN(UMAT,NBASIS,ECORE,ARR,BRR,G1)
         use constants, only: dp
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         USE UMatCache , only : UMatInd,UMAT2D,TUMAT2D
         use OneEInts, only: TMatind,TMat2D,TMATSYM,TSTARSTORE
         IMPLICIT NONE
         INTEGER NBASIS
         REAL*8 ECORE,ARR(NBASIS,2)
         REAL*8 UMAT(*)
         REAL*8 Z
         INTEGER*8 MASK,IND
         INTEGER I,J,K,L,BRR(NBASIS),G1(*),X,Y
         INTEGER NORB,NELEC,MS2,ORBSYM(1000),ISYM,SYMMAX
         INTEGER Counter(1:8),Index(1000)
         LOGICAL LWRITE,UHF
         NAMELIST /FCI/ NORB,NELEC,MS2,ORBSYM,ISYM,UHF
         LWRITE=.FALSE.
         UHF=.FALSE.
         OPEN(8,FILE='FCISYM',STATUS='OLD',FORM='FORMATTED')
         READ(8,FCI)
         CLOSE(8)
         OPEN(8,FILE='FCIDUMP',STATUS='OLD',FORM='UNFORMATTED')


!!Calculate a key, which will reorder the orbitals in symmetry order.
!!Counter is the cumulative number of symmetries, just to indicate the allowed current index of an orbital
!         Counter(:)=0
!         SymMax=1
!         do i=1,NORB
!             IF(ORBSYM(i).gt.SymMax) SymMax=ORBSYM(i)
!             Counter(ORBSYM(i))=Counter(ORBSYM(i))+1
!         enddo
!         do i=SymMax,2,-1
!             Counter(i)=1
!             do j=1,i-1
!                 Counter(i)=Counter(i)+Counter(j)
!             enddo
!         enddo
!         Counter(1)=1
!!Check that all orbitals are accounted for.
!         j=0
!         do i=1,NORB
!             IF(ORBSYM(i).eq.SymMax) j=j+1
!         enddo
!         IF((Counter(SymMax)+j-1).ne.NORB) THEN
!             CALL Stop_All("GETFCIBASIS","Counter incorrectly set up")
!         ENDIF
!!Now work out what the new indices are to be called
!         Index(1:NORB)=0
!         do i=1,norb
!             Index(i)=Counter(OrbSym(i))
!             Counter(OrbSym(i))=Counter(OrbSym(i))+1
!         enddo
!!Now need to sort OrbSym so that the symmetries correspond to the correct orbitals.
!         CALL NECI_SORTI(norb,OrbSym(1:norb))

         MASK=(2**16)-1
         !IND contains all the indices in an integer*8 - use mask of 16bit to extract them
101      READ(8,END=199) Z,IND
         L=iand(IND,MASK)
         IND=Ishft(IND,-16)
         K=iand(IND,MASK)
         IND=Ishft(IND,-16)
         J=iand(IND,MASK)
         IND=Ishft(IND,-16)
         I=iand(IND,MASK)
!         I=Index(I)
!         J=Index(J)
!         K=Index(K)
!         L=Index(L)
         
C.. Each orbital in the file corresponds to alpha and beta spinorbitalsa
         IF(I.EQ.0) THEN
C.. Core energy
            ECORE=Z
         ELSEIF(J.EQ.0) THEN
!C.. HF Eigenvalues
!            ARR(I*2-1,2)=Z
!            ARR(I*2,2)=Z
!            ARR(BRR(I*2-1),1)=Z
!            ARR(BRR(I*2),1)=Z
!            LWRITE=.TRUE.
         ELSEIF(K.EQ.0) THEN
C.. 1-e integrals
            IF(TSTARSTORE) THEN
! If TSTARSTORE, the one-el integrals are stored in symmetry classes, as spatial orbitals
                TMATSYM(TMatInd(2*J,2*I))=Z
            ELSE
C.. These are stored as spinorbitals (with elements between different spins being 0
                TMAT2D(2*I-1,2*J-1)=Z
                TMAT2D(2*I,2*J)=Z
                TMAT2D(2*J-1,2*I-1)=Z
                TMAT2D(2*J,2*I)=Z
            ENDIF
         ELSE
C.. 2-e integrals
C.. UMAT is stored as just spatial orbitals (not spinorbitals)
C..  we're reading in (IJ|KL), but we store <..|..> which is <IK|JL>
            IF(TSTARSTORE.and.TUMAT2D) THEN
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
                ELSE
                    UMAT(UMatInd(I,K,J,L,0,0))=Z
                ENDIF
            ELSEIF(TSTARSTORE.and.(.not.TUMAT2D)) THEN
                STOP 'Need UMAT2D with TSTARSTORE'
            ELSE
                UMAT(UMatInd(I,K,J,L,0,0))=Z
            ENDIF
!            UMAT(I,K,J,L)=Z
!            UMAT(J,L,I,K)=Z
!            UMAT(K,I,L,J)=Z
!            UMAT(L,J,K,I)=Z
!            UMAT(J,K,I,L)=Z
!            UMAT(I,L,J,K)=Z
!            UMAT(L,I,K,J)=Z
!            UMAT(K,J,L,I)=Z
         ENDIF
!         WRITE(14,'(1X,F20.12,4I3)') Z,I,J,K,L
         IF(I.NE.0) GOTO 101
199      CONTINUE
         CLOSE(8)
C.. If we've changed the eigenvalues, we write out the basis again
!         IF(LWRITE) THEN
!            WRITE(6,*) "1-electron energies have been read in."
!            CALL WRITEBASIS(6,G1,NBASIS,ARR,BRR)
!         ENDIF
         RETURN
      END
