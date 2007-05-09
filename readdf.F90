! lenrec is the number of auxiliary basis functions
SUBROUTINE InitDFBasis(nEl,nBasisMax,Len,lMs)
         use record_handler
         use HElement
         use UMatCache
         implicit none
         integer nEl,nBasisMax(5,3),Len,lMs
         parameter C_file='SAV_DFaSOL'
         parameter nolabel='        '
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
         use HElement
         use UMatCache
         implicit none
         parameter C_file='SAV_DFaSOL'
         parameter I_file='SAV_Ta_INT'
         parameter nolabel='        '
         character(3) file_status
         integer info,i
         integer nBasis,nOrbUsed,ierr
         real(dp) array(1000)

         file_status= 'ADD'
!.. We've already got C_file open
         call init_record_handler(I_file,file_status,info,printinfo=.TRUE.)
!.. lenrec is the number of auxiliary basis functions
!.. nrec is the number of pairs of orbitals.
         WRITE(6,*) "Basis Size:", nBasis/2
         WRITE(6,*) "Auxiliary Basis Size:", nAuxBasis

         tDFInts=.TRUE.
         Allocate(DFCoeffs(nAuxBasis,nBasisPairs),STAT=ierr)
         call MemAlloc(ierr,DFCoeffs,nBasisPairs*nAuxBasis,"DFCoeffs")
         Allocate(DFInts(nAuxBasis,nBasisPairs),STAT=ierr)
         call MemAlloc(ierr,DFInts,nBasisPairs*nAuxBasis,"DFInts")
         do i=1,nBasisPairs
            call read_record(C_file,i,nolabel,DFCoeffs(:,i),info)
            call read_record(I_file,i,nolabel,DFInts(:,i),info)
         enddo
         call leave_record_handler(C_file,info)
         call leave_record_handler(I_file,info)
         IF(nBasis.NE.nOrbUsed) THEN
! We allocate a small preliminary cache before freezing.
            call SetupUMatCache(nOrbUsed/2,.TRUE.)
         ELSE
            call SetupUMatCache(nOrbUsed/2,.FALSE.)
         ENDIF
         call SetupUMat2D_df
      END
      
!.. Get a 2-el integral.  a,b,c,d are indices.
      SUBROUTINE GetDF2EInt(a,b,c,d,res)
         use HElement 
         use UMatCache
         implicit none
         integer a,b,c,d
         integer i,GetDFIndex
         real*8 res
         res=0.D0
         do i=1,nAuxBasis
           res=res+DFCoeffs(i,GetDFIndex(a,c))*DFInts(i,GetDFIndex(b,d)) 
         enddo
!         WRITE(79,"(4I4,G)") a,b,c,d,res
!         WRITE(6,*) "D",a,b,c,d,res
      END
     
!.. return a DF pair index - i<j (although the pairs are ordered 11 21 22 31 32 33 41 42 ...
      INTEGER FUNCTION GetDFIndex(i,j)
         IMPLICIT NONE
         INTEGER I,J
         GetDFIndex=i+j*(j-1)/2
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
                  else
                     call GetDF2EInt(i,j,i,j,UMat2D(j,i))
                     call GetDF2EInt(i,j,j,i,UMat2D(i,j))
                  endif
                  WRITE(78,"(3I5,G)") i,j,0,UMat2D(i,j)
                  IF(i.ne.j) WRITE(78,"(3I5,G)") j,i,0,UMat2D(j,i)
               enddo
            enddo
            close(78)
         endif
      END
      SUBROUTINE ReadDalton1EIntegrals(G1,nBasis,TMat,Arr,Brr,ECore)
         implicit none
         include 'basis.inc'
         integer nBasis,Brr(nBasis),i,j
         real*8 Arr(nBasis,2),val,ECore,TMat(nBasis,nBasis)
         type(BasisFN) G1(nBasis)
         open(11,file='HONEEL',status='unknown')
         i=1
         call azZero(TMat,nBasis*nBasis)
         call iazZero(G1,nBasis*BasisFNSize)
         do while(i.ne.0)
            read(11,*) i,j,val
!"(2I5,F)"
            if(i.eq.0) then
               ECore=val
            elseif(j.ne.0) then
               TMat(i*2-1,j*2-1)=val
               TMat(i*2,j*2)=val
               TMat(j*2-1,i*2-1)=val
               TMat(j*2,i*2)=val
            endif
         enddo
         close(11)
         do i=1,nBasis
            G1(i)%Ms=1-2*iand(i,1)
            G1(i)%Sym%s=1
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
         real*8 Arr(nBasis,2),val,ECore,TMat(nBasis,nBasis)
         type(BasisFN) G1(nBasis)
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
            G1(i)%Sym%s=1
            Brr(i)=i
         enddo
      END 
