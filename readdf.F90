! lenrec is the number of auxiliary basis functions
SUBROUTINE InitDFBasis(nEl,nBasisMax,Len,lMs)
         use record_handler
         use HElement
         implicit none
         include 'umatcache.inc'
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
         implicit none
         INCLUDE 'umatcache.inc'
         parameter C_file='SAV_DFaSOL'
         parameter I_file='SAV_Ta_INT'
         parameter nolabel='        '
         character(3) file_status
         integer info,i
         integer nBasis,nOrbUsed
         real(dp) array(1000)

         file_status= 'ADD'
!.. We've already got C_file open
         call init_record_handler(I_file,file_status,info,printinfo=.TRUE.)
!.. lenrec is the number of auxiliary basis functions
!.. nrec is the number of pairs of orbitals.
         WRITE(6,*) "Basis Size:", nBasis/2
         WRITE(6,*) "Auxiliary Basis Size:", nAuxBasis

         call Memory(IP_DFCoeffs,nBasisPairs*nAuxBasis,"DFCoeffs")
         call Memory(IP_DFInts,nBasisPairs*nAuxBasis,"DFInts")
         do i=1,nBasisPairs
            call read_record(C_file,i,nolabel,DFCoeffs(:,i),info)
            call read_record(I_file,i,nolabel,DFInts(:,i),info)
         enddo
         call leave_record_handler(C_file,info)
         call leave_record_handler(I_file,info)
         call SetupUMatCache(nOrbUsed/2)
      END
      
!.. Get a 2-el integral.  a,b,c,d are indices.
      SUBROUTINE GetDF2EInt(a,b,c,d,res)
         use HElement 
         implicit none
         include 'umatcache.inc'
         integer a,b,c,d
         integer i,GetDFIndex
         real*8 res
         res=0.D0
         do i=1,nAuxBasis
           res=res+DFCoeffs(i,GetDFIndex(a,c))*DFInts(i,GetDFIndex(b,d)) 
         enddo
      END
     
!.. return a DF pair index - i<j (although the pairs are ordered 11 21 22 31 32 33 41 42 ...
      INTEGER FUNCTION GetDFIndex(i,j)
         IMPLICIT NONE
         INTEGER I,J
         GetDFIndex=i+j*(j-1)/2
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
            if(i.eq.0) then
               ECore=val
            else
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
            Arr(i,1)=TMat(i,i)
            Arr(i,2)=TMat(i,i)
            Brr(i)=i
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
            if(i.eq.j.and.i.ne.0) then
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
