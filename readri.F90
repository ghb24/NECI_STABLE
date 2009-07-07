subroutine InitRIBasis(nEl,nBasisMax,Len,lMs)
   USE HElem
! lenrec is the number of auxiliary basis functions
   use UMatCache
   implicit none
   integer nEl,nBasisMax(5,*),Len,lMs
   integer info,lenrec,nrec,i
   integer nBasis
   integer*8 nAb,nB
   integer record_length
   OPEN(29,file='RIINTDUMP',status='old',FORM='UNFORMATTED',access='DIRECT',recl=record_length(8))
!.. The first element is the number of aux basis fns.
!.. The second element is the number of basisfunctions.
   READ(29,rec=1) nAb
   READ(29,rec=2) nB
   write(6,*) nAb,nB
   nAuxBasis=nAb
   nBasis=nB
   WRITE(6,*) "Q-Chem auxiliary basis", nAuxBasis, " basis functions:", nBasis
   nBasisMax(1:5,1:3)=0
   Len=2*nBasis 
!.. Note that it's a read in basis.
   nBasisMax(3,3)=1
   nBasisMax(4,1)=-1
   nBasisMax(4,2)=1
!.. Correspond to ISS=0
   nBasisMax(1,3)=0
!.. Setup Max Sym 
   nBasisMax(5,2)=0
   CLOSE(29)
END

SUBROUTINE GetRI2EInt(a,b,c,d,res)
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
! DFOVERLAP        1 - (ij|u|ab)= (ij|u|P)(P|ab)
   do i=1,nAuxBasis
     res=res+DFCoeffs(i,x)*DFInts(i,y)
   enddo
end

SUBROUTINE ReadRI2EIntegrals(nBasis,nOrbUsed)
   use UMatCache
   IMPLICIT NONE
   INTEGER nBasis,nOrbUsed
   IF(nBasis.NE.nOrbUsed) THEN
! We allocate a small preliminary cache before freezing.
      WRITE(6,*) "Setting up pre-freezing UMatCache"
      call SetupUMatCache(nOrbUsed/2,.TRUE.)
   ELSE
      call SetupUMatCache(nOrbUsed/2,.FALSE.)
   ENDIF
   call SetupUMat2D_df()
!Integrals are actually read in READFCIINT.  later on
END

SUBROUTINE ReadRIIntegrals(nBasis,nOrbUsed)
   use UMatCache
   use global_utilities
   IMPLICIT NONE
   character(*), parameter :: t_r='ReadRIIntegrals'
   INTEGER nBasis,nOrbUsed
   integer i,j,k,ilast,onints,nints,ierr,Q
   real*8 val
   integer*8 nA,nB
   integer GetDFIndex
   integer record_length
   WRITE(6,*) "Reading QChem C Matrices"
   OPEN(29,file='RIINTDUMP',status='old',FORM='UNFORMATTED',access='DIRECT',recl=record_length(8))
   read(29,rec=1) nA
   read(29,rec=2) nB
   onints=2
   nints=(nBasis/2*(nBasis/2+1))/2  !/2 to convert to states
   Allocate(DFInts(nAuxBasis,nints),STAT=ierr)
   call LogMemAlloc("DFInts",nints*nAuxBasis,8,t_r,tagDFInts,ierr)
   WRITE(6,*) nints,nAuxBasis
   DO I=1,nBasis/2
      do Q=1,nAuxBasis
         do J=1,nBasis/2
            onints=onints+1
            read(29,rec=onints) val
            DFInts(Q,GetDFIndex(j,i))=val
         enddo
      enddo
   enddo
   WRITE(6,*) "Over"
   close(29)
   iDFMethod=3
   CALL ReadRI2EIntegrals(nBasis,nOrbUsed)
   WRITE(6,*) "DFM:",iDFMethod
end
