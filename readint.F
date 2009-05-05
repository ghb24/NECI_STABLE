      SUBROUTINE INITFROMFCID(NEL,NBASISMAX,LEN,LMS,TBIN)
         IMPLICIT NONE
         INTEGER NEL,nBasisMax(5,*),LEN,LMS
         INTEGER NORB,NELEC,MS2,ORBSYM(1000),ISYM
         LOGICAL TBIN,exists,UHF
         CHARACTER*3 :: fmat
         NAMELIST /FCI/ NORB,NELEC,MS2,ORBSYM,ISYM,UHF
         UHF=.FALSE.
         fmat='NO'
         IF(TBIN) THEN
            INQUIRE(FILE='FCISYM',EXIST=exists)
            IF(.not.exists) STOP 'FCISYM file does not exist'
            INQUIRE(FILE='FCIDUMP',EXIST=exists,FORMATTED=fmat)
            IF(.not.exists) STOP 'FCIDUMP file does not exist'
!            IF(fmat=='YES') THEN
!                STOP 'FCIDUMP is not unformatted, but TBIN true'
!            ENDIF
            OPEN(8,FILE='FCISYM',STATUS='OLD',FORM='FORMATTED')
            READ(8,FCI)
         ELSE
            INQUIRE(FILE='FCIDUMP',EXIST=exists,UNFORMATTED=fmat)
            IF(.not.exists) STOP 'FCIDUMP file does not exist'
!            IF(fmat=='YES') THEN
!                STOP 'FCIDUMP is not formatted, but TBIN false'
!            ENDIF
            OPEN(8,FILE='FCIDUMP',STATUS='OLD',FORM='FORMATTED')
            READ(8,FCI)
         ENDIF
         CLOSE(8)
         IF(NELEC.NE.NEL)  WRITE(6,*) 
     &      '*** WARNING: NEL in FCIDUMP differs from input file ***'
C         NEL=NELEC
         IF(LMS.NE.MS2)   WRITE(6,*)
     &      '*** WARNING: LMS in FCIDUMP differs from input file ***'
         LMS=MS2
         WRITE(6,*) ' NUMBER OF ELECTRONS : ' , NEL
         WRITE(6,*) ' BASIS MS : ' , LMS
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
         IF(UHF) then
            NBASISMAX(2,3)=1
         else
            NBASISMAX(2,3)=2
         endif

C.. SHow that there's no momentum conservation
         NBASISMAX(3,3)=1
         RETURN
      END
      SUBROUTINE GETFCIBASIS(NBASISMAX,ARR,BRR,G1,LEN,TBIN)
         use SystemData, only: BasisFN,BasisFNSize,Symmetry
         IMPLICIT NONE
         INTEGER nBasisMax(5,*),BRR(LEN),LEN
         TYPE(BasisFN) G1(LEN)
         REAL*8 ARR(LEN,2)
         REAL*8 Z
         INTEGER*8 IND,MASK
         INTEGER I,J,K,L
         INTEGER ISYMNUM,ISNMAX,SYMMAX
         INTEGER NORB,NELEC,MS2,ORBSYM(1000),ISYM,ISPINS,ISPN
         INTEGER Counter(1:8),Index(1000)
         LOGICAL TBIN,UHF
         NAMELIST /FCI/ NORB,NELEC,MS2,ORBSYM,ISYM,UHF
         UHF=.FALSE.
         IF(TBIN) THEN
            OPEN(8,FILE='FCISYM',STATUS='OLD',FORM='FORMATTED')
            READ(8,FCI)
            CLOSE(8)
            OPEN(8,FILE='FCIDUMP',STATUS='OLD',FORM='UNFORMATTED')
         ELSE
            OPEN(8,FILE='FCIDUMP',STATUS='OLD',FORM='FORMATTED')
            READ(8,FCI)
         ENDIF
         ISYMNUM=0
         ISNMAX=0
         ISPINS=2
         IF(UHF) ISPINS=1
        
         IF(LEN.NE.ISPINS*NORB) STOP 'LEN .NE. NORB in GETFCIBASIS'
         G1(1:LEN)=BasisFN((/0,0,0/),0,Symmetry(0))
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



         IF(TBIN) THEN
            IF(UHF) STOP 'UHF Bin read not functional'
            MASK=(2**16)-1
            
            !IND contains all the indices in an integer*8 - use mask of 16bit to extract them
2           READ(8,END=99) Z,IND
            L=iand(IND,MASK)
            IND=Ishft(IND,-16)
            K=iand(IND,MASK)
            IND=Ishft(IND,-16)
            J=iand(IND,MASK)
            IND=Ishft(IND,-16)
            I=iand(IND,MASK)

!            I=Index(I)
!            J=Index(J)
!            K=Index(K)
!            L=Index(L)
            
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
         ELSE
1           READ(8,'(1X,G20.12,4I3)',END=99) Z,I,J,K,L
!            I=Index(I)
!            J=Index(J)
!            K=Index(K)
!            L=Index(L)
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
                DO ISPN=1,ISPINS
                    ARR(ISPINS*I-ISPN+1,1)=Z
                ENDDO
            ELSEIF(I.NE.0.AND.K.EQ.0.AND.J.EQ.0) THEN
                DO ISPN=1,ISPINS
                    ARR(ISPINS*I-ISPN+1,1)=Z
                ENDDO
            ENDIF
C.. At the moment we're ignoring the core energy
            IF(I.NE.0) GOTO 1
         ENDIF
99       CONTINUE
         SYMMAX=1
         DO I=1,NORB
            DO ISPN=1,ISPINS
                BRR(ISPINS*I-ISPN+1)=ISPINS*I-ISPN+1
                G1(ISPINS*I-ISPN+1)%Sym%s=ORBSYM(I)-1
C.. set momentum to 0
                G1(ISPINS*I-ISPN+1)%k(1)=0
                G1(ISPINS*I-ISPN+1)%k(2)=0
                G1(ISPINS*I-ISPN+1)%k(3)=0
                G1(ISPINS*I-ISPN+1)%Ms=-MOD(ISPINS*I-ISPN+1,2)*2+1
               IF(SYMMAX.lt.ORBSYM(I)) SYMMAX=ORBSYM(I)
            ENDDO
         ENDDO
         IF(SYMMAX.EQ.1) THEN
         ELSEIF(SYMMAX.LE.2) THEN
         ELSEIF(SYMMAX.LE.4) THEN
            SYMMAX=4
         ELSEIF(SYMMAX.LE.8) THEN
            SYMMAX=8
         ENDIF
        WRITE(6,"(A,I3)") "Maximum number of symmetries: ",SYMMAX
         NBASISMAX(1,1)=0
         NBASISMAX(1,2)=0
         NBASISMAX(5,1)=0
         NBASISMAX(5,2)=SYMMAX-1
         NBASISMAX(2,1)=0
         NBASISMAX(2,2)=0
         CLOSE(8)
         RETURN
      END
      SUBROUTINE READFCIINT(UMAT,NBASIS,ECORE,ARR,BRR,G1)
         USE HElem
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SystemData, only: UMatEps,tUMatEps
         USE UMatCache, only: UMatInd,UMAT2D,TUMAT2D
         use OneEInts, only: TMatind,TMat2D,TMATSYM,TSTARSTORE
         IMPLICIT NONE
         INTEGER NBASIS,ZeroedInt,NonZeroInt
         REAL*8 ECORE,ARR(NBASIS,2)
         REAL*8 UMAT(*)
         REAL*8 Z
         INTEGER I,J,K,L,BRR(NBASIS),G1(*),X,Y
         INTEGER NORB,NELEC,MS2,ORBSYM(1000),ISYM,SYMMAX
         LOGICAL LWRITE,UHF
         INTEGER ISPINS,ISPN,ISPN2
         INTEGER Counter(1:8),Index(1000)
         NAMELIST /FCI/ NORB,NELEC,MS2,ORBSYM,ISYM,UHF
         LWRITE=.FALSE.
         UHF=.FALSE.
         ZeroedInt=0
         NonZeroInt=0
         
         OPEN(8,FILE='FCIDUMP',STATUS='OLD')
         READ(8,FCI)
         ISPINS=2
         IF(UHF) ISPINS=1

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



101      READ(8,'(1X,G20.12,4I3)',END=199) Z,I,J,K,L
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
                DO ISPN=1,ISPINS
                   TMAT2D(ISPINS*I-ISPN+1,ISPINS*J-ISPN+1)=Z
                   TMAT2D(ISPINS*J-ISPN+1,ISPINS*I-ISPN+1)=Z
                enddo
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
199      CONTINUE
         IF(tUMatEps) THEN
!Write out statistics from zeroed integrals.
             WRITE(6,"(A,G20.10,A)") "*** Zeroing all two-electron "    &
     &          //"integrals with a magnitude of over ",UMatEps," ***" 
             WRITE(6,*) ZeroedInt+NonZeroInt," 2E integrals read in..."
             WRITE(6,*) ZeroedInt," integrals zeroed..."
             WRITE(6,*) REAL(100*ZeroedInt)/REAL(NonZeroInt+ZeroedInt), &
     %          " percent of 2E integrals zeroed."
         ENDIF
         CLOSE(8)
C.. If we've changed the eigenvalues, we write out the basis again
!         IF(LWRITE) THEN
!            WRITE(6,*) "1-electron energies have been read in."
!            CALL WRITEBASIS(6,G1,NBASIS,ARR,BRR)
!         ENDIF
         RETURN
      END

      !This is a copy of the routine above, but now for reading in binary files of integrals
      SUBROUTINE READFCIINTBIN(UMAT,NBASIS,ECORE,ARR,BRR,G1)
         USE HElem
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
