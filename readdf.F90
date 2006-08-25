! lenrec is the number of auxiliary basis functions
SUBROUTINE InitDFBasis(nEl,nBasisMax,Len,lMs)
         use record_handler
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
      END
!  A file to read in density fitted two-electron integrals from SITUS/Dalton.
!.. Also will read in one-electron integrals
!.. Requires various modules from SITUS for the file reading
      SUBROUTINE ReadDF2EIntegrals(nBasis)
         use precision
         use record_handler
         implicit none
         INCLUDE 'umatcache.inc'
         parameter C_file='SAV_DFaSOL'
         parameter I_file='SAV_Ta_INT'
         parameter nolabel='        '
         character(3) file_status
         integer info,lenrec,nrec,i
         integer nBasis
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
         do i=1,nrec
            call read_record(C_file,i,nolabel,DFCoeffs(:,i),info)
            call read_record(I_file,i,nolabel,DFInts(:,i),info)
         enddo
         call leave_record_handler(C_file,info)
         call leave_record_handler(I_file,info)
      END
      
!.. Get a 2-el integral.  a and b are pair indices.
      SUBROUTINE GetDF2EInt(a,b,res)
         implicit none
         include 'umatcache.inc'
         integer a,b
         integer i
         real*8 res
         res=0.D0
         do i=1,nAuxBasis
           res=res+DFCoeffs((a-1)*nAuxBasis)*DFInts(i,b) 
         enddo
      END
      
      SUBROUTINE ReadDalton1EIntegrals
      END
   
      SUBROUTINE InitDaltonBasis
      END
