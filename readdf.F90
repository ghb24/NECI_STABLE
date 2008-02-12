! lenrec is the number of auxiliary basis functions
SUBROUTINE InitDFBasis(nEl,nBasisMax,Len,lMs)
         use record_handler
         USE HElem
         use UMatCache
         implicit none
         integer nEl,nBasisMax(5,3),Len,lMs
         character(*), parameter :: C_file='SAV_D____a'
         character(*), parameter :: nolabel='        '
         character(3) file_status
         integer info,lenrec,nrec,i
         integer nBasis

         file_status= 'ADD'
         call init_record_handler(C_file,file_status,info)
         call query_record_handler(C_file,info,file_status,lenrec,nrec,.TRUE.)
!.. lenrec is the number of auxiliary basis functions
!.. nrec is the number of pairs of orbitals.
         nBasisPairs=nrec
         nBasis=int(dsqrt(nBasisPairs*2.D0))
         nAuxBasis=lenrec
         WRITE(6,*) "DALTON/SITUS basis.", nBasis, " basis functions."
         call iAzZero(nBasisMax,15)
         Len=2*nBasis 
!.. Note that it's a read in basis.
         nBasisMax(3,3)=1
         nBasisMax(4,1)=-1
         nBasisMax(4,2)=1
!.. Correspond to ISS=0
         nBasisMax(1,3)=0
!.. Setup Max Sym 
         nBasisMax(5,2)=0
      END
!  A file to read in density fitted two-electron integrals from SITUS/Dalton.
!.. Also will read in one-electron integrals
!.. Requires various modules from SITUS for the file reading
      SUBROUTINE ReadDF2EIntegrals(nBasis,nOrbUsed)
         use precision
         use record_handler
         USE HElem
         use UMatCache
         implicit none
         character(*), parameter :: C_file='SAV_D____a'
         character(*), parameter :: I_file='SAV_T____a'
         character(*), parameter :: S_file='SAV_S____a'
         character(*), parameter :: nolabel='        '
         character(3) file_status
         integer info,i,j,k
         real*8 r
         integer nBasis,nOrbUsed,ierr
         real(dp) array(1000)

         WRITE(6,*) "Opening Density fitting matrix files"
         file_status= 'ADD'
!.. We've already got C_file open
         call init_record_handler(I_file,file_status,info,printinfo=.TRUE.)
         call init_record_handler(S_file,file_status,info,printinfo=.TRUE.)
!.. lenrec is the number of auxiliary basis functions
!.. nrec is the number of pairs of orbitals.
         WRITE(6,*) "Basis Size:", nBasis/2
         WRITE(6,*) "Auxiliary Basis Size:", nAuxBasis

         tDFInts=.TRUE.
         Allocate(DFCoeffs(nAuxBasis,nBasisPairs),STAT=ierr)
         call MemAlloc(ierr,DFCoeffs,nBasisPairs*nAuxBasis,"DFCoeffs")
         Allocate(DFInts(nAuxBasis,nBasisPairs),STAT=ierr)
         call MemAlloc(ierr,DFInts,nBasisPairs*nAuxBasis,"DFInts")
         Allocate(DFFitInts(nAuxBasis,nAuxBasis),STAT=ierr)
         call MemAlloc(ierr,DFFitInts,nAuxBasis*nAuxBasis,"DFFitInts")
         do i=1,nBasisPairs
            call read_record(C_file,i,nolabel,DFCoeffs(:,i),info)
            call read_record(I_file,i,nolabel,DFInts(:,i),info)
         enddo
         do i=1,nAuxBasis
            call read_record(S_file,i,nolabel,DFFitInts(:,i),info)
         enddo
         call leave_record_handler(C_file,info)
         call leave_record_handler(I_file,info)
         call leave_record_handler(S_file,info)

         select case(iDFMethod)
         case(0)
            STOP "No DF, but ReadDF2EIntegrals called."
         case(1)
            WRITE(6,*) "DFMETHOD DFOVERLAP        1 - (ij|u|ab)= (ij|u|P)(P|ab)"
         case(2)
            WRITE(6,*) "DFMETHOD DFOVERLAP2NDORD  2 - (ij|u|ab)= (ij|u|P)(P|ab)+(ij|P)(P|u|ab)-(ij|P)(P|u|Q)(Q|ab)"
         case(3)
            WRITE(6,*) "DFMETHOD DFOVERLAP2       3 - (ij|u|ab)= (ij|P)(P|u|Q)(Q|ab)"
         case(4)
            WRITE(6,*) "DFMETHOD DFCOULOMB        4 - (ij|u|ab)= (ij|u|P)[(P|u|Q)^-1](Q|u|ab)"
         case default
            WRITE(6,*) "Unknown DF Method: ",iDFMethod
            STOP "Unknown DF Method"
         end select

         if(iDFMethod.gt.2) then
            Allocate(DFInvFitInts(nAuxBasis,nAuxBasis),STAT=ierr)
            call MemAlloc(ierr,DFInvFitInts,nAuxBasis*nAuxBasis,"DFInvFitInts")
         endif
         select case(iDFMethod)
         case(3)
            call DFCalcInvFitInts(0.5D0)
            WRITE(6,*) "Calculating B matrix" 
            do i=1,nBasisPairs
               do j=1,nAuxBasis
                  r=0
                  do k=1,nAuxBasis
                     r=r+DFInvFitInts(j,k)*DFCoeffs(k,i)
                  enddo
                  DFInts(j,i)=r
               enddo
            enddo
!DFInts now contains B_ij,P = sum_Q (ij|Q)[(Q|u|P)^1/2]
         case(4)
            call DFCalcInvFitInts(-0.5D0)
            WRITE(6,*) "Calculating B matrix" 
            do i=1,nBasisPairs
               do j=1,nAuxBasis
                  r=0
                  do k=1,nAuxBasis
                     r=r+DFInvFitInts(j,k)*DFInts(k,i)
                  enddo
                  DFCoeffs(j,i)=r
               enddo
            enddo
!DFCoeffs now contains B_ij,P = sum_Q (ij|u|Q)[(Q|u|P)^-1/2]
         end select
         IF(nBasis.NE.nOrbUsed) THEN
! We allocate a small preliminary cache before freezing.
            call SetupUMatCache(nOrbUsed/2,.TRUE.)
         ELSE
            call SetupUMatCache(nOrbUsed/2,.FALSE.)
         ENDIF
         call SetupUMat2D_df
      END
      
!.. Get a 2-el integral.  a,b,c,d are indices. <ab|1/r12|cd>
!DFCoeffs(x,yz) is (x|yz)
!DFInts(x,yz) is (x|u|yz)
!DFFitInts(x,y) is (x|u|y)
      SUBROUTINE GetDF2EInt(a,b,c,d,res)
         USE HElem 
         use UMatCache
         implicit none
         integer a,b,c,d
         integer i,j,GetDFIndex
         integer x,y
         real*8 res
         res=0.D0
         x=GetDFIndex(a,c)
         y=GetDFIndex(b,d)
!CC         write(6,"(6I4)") a,b,c,d,x,y
         select case(iDFMethod)
! 0 - no DF
         case(0)
            STOP "DF Method 0 specified - no DF, but DF called."
         case(1)
! DFOVERLAP        1 - (ij|u|ab)= (ij|u|P)(P|ab)
            do i=1,nAuxBasis
              res=res+DFCoeffs(i,x)*DFInts(i,y)
            enddo
         case(2)
! DFOVERLAP2NDORD  2 - (ij|u|ab)= (ij|u|P)(P|ab)+(ij|P)(P|u|ab)-(ij|P)(P|u|Q)(Q|ab)
            do i=1,nAuxBasis
              res=res+DFCoeffs(i,x)*DFInts(i,y)+DFCoeffs(i,y)*DFInts(i,x)
              do j=1,nAuxBasis
                  res=res-DFCoeffs(i,x)*DFCoeffs(j,y)*DFFitInts(i,j)
               enddo
            enddo
         case(3)
! DFOVERLAP2       3 - (ij|u|ab)= (ij|P)(P|u|Q)(Q|ab)

!DFInts actually contains B_ij,P = sum_Q (ij|Q)[(Q|u|P)^1/2]
            do i=1,nAuxBasis
              res=res+DFInts(i,x)*DFInts(i,y)
            enddo
         case(4)
! DFCOULOMB        4 - (ij|u|ab)= (ij|u|P)[(P|u|Q)^-1](Q|u|ab)
!DFCoeffs actually contains B_ij,P = sum_Q (ij|u|Q)[(Q|u|P)^-1/2]
            do i=1,nAuxBasis
              res=res+DFCoeffs(i,x)*DFCoeffs(i,y)
            enddo
         end select
      END


     
!.. return a DF pair index - i<j (although the pairs are ordered 11 21 22 31 32 33 41 42 ...
      INTEGER FUNCTION GetDFIndex(i,j)
         IMPLICIT NONE
         INTEGER I,J
         if(i.lt.j) then
            GetDFIndex=i+j*(j-1)/2
         else
            GetDFIndex=j+i*(i-1)/2
         endif
      END 
      SUBROUTINE ReadDalton2EIntegrals(nBasis,UMat2D,tUMat2D)
         implicit none
         include 'basis.inc'
         integer nBasis,i,j,k,ilast
         real*8 val,UMat2D(nBasis,nBasis)
         logical tUMat2D
         tUMat2D=.false.
         open(11,file='HONEEL',status='unknown')
         i=1
         do while(i.ne.0)
            read(11,*) i,j,val
         enddo
         i=1
         ilast=0
         do while(i.ne.0.and.i.le.nBasis.and.j.le.nBasis)
            read(11,*,end=20) i,j,k,val
            if(i.ne.0.and.i.le.nBasis.and.j.le.nBasis) then
               UMat2D(i,j)=val
               ilast=i
            endif
         enddo
         tUMat2D=.true.
20       close(11)
         IF(tUMat2D) THEN
            WRITE(6,*) "Read in 2-index 2-electron integrals up to orbital ", ilast*2
         ELSE
            WRITE(6,*) "Have not read in 2-index 2-electron integrals."
         ENDIF
         IF(ilast.lt.nBasis) THEN
            WRITE(6,*) "Calculating remaining 2-index 2-electron integrals with density fitting."
            OPEN(78,file='UMAT',status='UNKNOWN')
            do i=ilast+1,nBasis
               do j=1,i
                  if(i.lt.j) then
                     call GetDF2EInt(i,j,i,j,UMat2D(i,j))
                     call GetDF2EInt(i,j,j,i,UMat2D(j,i))
                     WRITE(78,"(3I5,G)") i,j,0,UMat2D(i,j)
                     IF(i.ne.j) WRITE(78,"(3I5,G)") j,i,0,UMat2D(j,i)
                  else
                     call GetDF2EInt(i,j,i,j,UMat2D(j,i))
                     call GetDF2EInt(i,j,j,i,UMat2D(i,j))
                     WRITE(78,"(3I5,G)") i,j,0,UMat2D(j,i)
                     IF(i.ne.j) WRITE(78,"(3I5,G)") j,i,0,UMat2D(i,j)
                  endif
               enddo
            enddo
            close(78)
         endif
      END
      SUBROUTINE ReadDalton1EIntegrals(G1,nBasis,Arr,Brr,ECore)
         USE HElem
         USE UMatcache , only : TMATind,TMAT2D,TMATSYM,TSTARSTORE
         implicit none
         include 'basis.inc'
         integer nBasis,Brr(nBasis),i,j
         real*8 Arr(nBasis,2),val,ECore
         type(BasisFN) G1(nBasis)
         integer*8 TotSymRep
         open(11,file='HONEEL',status='unknown')
         i=1
         !call azZero(TMat,nBasis*nBasis)
         call iazZero(G1,nBasis*BasisFNSize)
         do while(i.ne.0)
            read(11,*) i,j,val
!"(2I5,F)"
            if(i.eq.0) then
               ECore=val
            elseif(j.ne.0) then
                IF(TSTARSTORE) THEN
                    TMatSYM(TMATInd(i*2-1,j*2-1))=HElement(val)
                    TMatSYM(TMATInd(i*2,j*2))=HElement(val)
                    TMatSYM(TMATInd(j*2-1,i*2-1))=HElement(val)
                    TMatSYM(TMATInd(j*2,i*2))=HElement(val)
                ELSE
                    TMat2D(i*2-1,j*2-1)=HElement(val)
                    TMat2D(i*2,j*2)=HElement(val)
                    TMat2D(j*2-1,i*2-1)=HElement(val)
                    TMat2D(j*2,i*2)=HElement(val)
                ENDIF
            endif
         enddo
         close(11)
         do i=1,nBasis
            G1(i)%Ms=1-2*iand(i,1)
            G1(i)%Sym%s=TotSymRep()
!  We've already read in and ordered the Energies
!            Arr(i,1)=TMat(i,i)
!            Arr(i,2)=TMat(i,i)
!            Brr(i)=i
         enddo
      END
      SUBROUTINE InitDaltonBasis(nBasisMax,Arr,Brr,G1,nBasis)
         implicit none
         include 'basis.inc'
         integer nBasis,Brr(nBasis),i,j,nBasisMax(5,3)
         real*8 Arr(nBasis,2),val,ECore
         type(BasisFN) G1(nBasis)
         integer*8 TotSymRep
         open(11,file='HONEEL',status='unknown')
         i=1
         call iazZero(G1,nBasis*BasisFNSize)
         do while(i.ne.0)
            read(11,*) i,j,val
            if(j.eq.0.and.i.ne.0) then
               Arr(i*2-1,1)=val
               Arr(i*2-1,2)=val
               Arr(i*2,1)=val
               Arr(i*2,2)=val
            endif
         enddo
         close(11)
         do i=1,nBasis
            G1(i)%Ms=1-2*iand(i,1)
            G1(i)%Sym%s=TotSymRep()
            Brr(i)=i
         enddo
      END 


!.. Get a 2-el integral.  a,b,c,d are indices. <ab|1/r12|cd>
!DFCoeffs(x,yz) is (x|yz)
!DFInts(x,yz) is (x|u|yz)
!DFFitInts(x,y) is (x|u|y)
!This is slower but calculates more accurately.
      SUBROUTINE GetDF2EInt2OrderOverlap(a,b,c,d,res)
         USE HElem 
         use UMatCache
         implicit none
         integer a,b,c,d
         integer i,GetDFIndex
         integer x,y,j
         real*8 res,res1,res2,res3
         res=0.D0
         x=GetDFIndex(a,c)
         y=GetDFIndex(b,d)
!  (ab|u|cd)=sum_PQ (ab|P)(P|u|Q)(Q|cd)
         res1=0
         res2=0
         res3=0
         do i=1,nAuxBasis
           do j=1,nAuxBasis
               res=res+DFCoeffs(i,x)*DFCoeffs(j,y)*DFFitInts(i,j)
            enddo
         enddo
!         WRITE(79,"(4I4,G)") a,b,c,d,res
!         WRITE(6,*) "D2",a,b,c,d,res,res1,res2,res3
      END

      SUBROUTINE GetDF2EInt2OrderCoulomb(a,b,c,d,res)
         USE HElem 
         use UMatCache
         implicit none
         integer a,b,c,d
         integer i,GetDFIndex
         integer x,y,j
         real*8 res,res1,res2,res3
         res=0.D0
         x=GetDFIndex(a,c)
         y=GetDFIndex(b,d)
!  (ab|u|cd)=sum_PQ (ab|u|P)[(P|u|Q)^-1](Q|u|cd)
         res1=0
         res2=0
         res3=0
         do i=1,nAuxBasis
           do j=1,nAuxBasis
               res=res+DFInts(i,x)*DFInts(j,y)*DFInvFitInts(i,j)
            enddo
         enddo
!         WRITE(79,"(4I4,G)") a,b,c,d,res
!         WRITE(6,*) "D2",a,b,c,d,res,res1,res2,res3
      END
      SUBROUTINE DFCalcInvFitInts(dPower)
         USE HElem
         use UMatCache
         implicit none
         Real*8,Pointer :: M(:,:) !(nAuxBasis,nAuxBasis)
         Real*8 Eigenvalues(nAuxBasis),r,dPower
         Real*8 Work(3*nAuxBasis)
         integer Workl,info
         Integer iSub
         Integer i,j,ierr,k,iMinEigv
         CALL TISET('DFInvFitIn',ISUB)
         Allocate(M(nAuxBasis,nAuxBasis),STAT=ierr)
         call MemAlloc(ierr,M,nAuxBasis*nAuxBasis,"M-DFInvFitInts")
         call azzero(M,nAuxBasis*nAuxBasis)
         do i=1,nAuxBasis
            do j=1,i
               M(i,j)=DFFitInts(i,j)
            enddo
         enddo
         Workl=3*nAuxBasis
         WRITE(6,*) "Diagonalizing (P|u|Q)"
         CALL DSYEV('V','L',nAuxBasis,M,nAuxBasis,Eigenvalues,WORK,WORKL,INFO)
         IF(INFO.NE.0) THEN
            WRITE(6,*) 'DYSEV error: ',INFO
            STOP
         ENDIF
         iMinEigv=1
         if(dPower.lt.0) then
            do i=1,nAuxBasis
               if(Eigenvalues(i).lt.1d-10) iMinEigv=i+1
            enddo
            write(6,*) "Ignoring ",iMinEigv-1," eigenvalues <1d-10"
         endif
!M now contains eigenvectors.  Eigenvector i is in M(:,I)
! A=U^T L U (L is the matrix of eigenvalues on the diagonal)
! A^-1 = U^T L^-1 U
         WRITE(6,*) "Calculating (P|u|Q)^",dPower
         do i=1,nAuxBasis
            do j=1,i
               r=0
               do k=iMinEigv,nAuxBasis
                  r=r+(Eigenvalues(k)**dPower)*M(j,k)*M(i,k)
               end do
               DFInvFitInts(i,j)=r
               DFInvFitInts(j,i)=r
            enddo
         enddo
         call MemDealloc(M)         
         Deallocate(M)
         CALL TIHALT('DFInvFitIn',ISUB)
      END
