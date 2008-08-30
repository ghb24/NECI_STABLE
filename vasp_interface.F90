module vasp_interface

implicit none
public
save
integer, parameter :: q=kind(0.d0)

integer, allocatable :: KPntInd(:)
real(q), allocatable :: eigv(:),kpnts(:,:)
real(q) :: xi

integer :: nKP,nStates,KPMsh(3)

contains

   subroutine NECIReceiveVASPData(vasp_nbands,vasp_nkpts,vasp_kpnts,vasp_kpmsh,vasp_xi,vasp_kpntind,vasp_eigv,vasp_umat2D,vasp_nEl)
      use HElem
      use SystemData, only: nEl
      use UMatCache, only: SetupUMatCache,UMat2D,tagUMat2D,tUMat2D
      implicit none
      integer :: vasp_nbands,vasp_nkpts,vasp_kpntind(vasp_nbands*vasp_nkpts),vasp_kpmsh(3)
      real(q) :: vasp_kpnts(3,vasp_nkpts),vasp_xi,vasp_eigv(vasp_nbands*vasp_nkpts),vasp_nEl
      complex(q) :: vasp_umat2d(vasp_nbands*vasp_nkpts,vasp_nbands*vasp_nkpts)
      integer :: i,j,ik,ierr

      nKP=vasp_nkpts
      nStates=vasp_nbands*vasp_nkpts

      write (6,*) 'nkp',nKP
      write (6,*) 'kp'
      do i=1,nKP
         write (6,*) vasp_kpnts(:,i)
      end do 
      nEl=nint(vasp_nEl*nKP)

      allocate(KPNTInd(nStates))
      allocate(eigv(nStates))
      allocate(kpnts(3,nKP))

      KPMsh(:)=vasp_kpmsh(:)
      KPNTInd(:)=vasp_kpntind(:)
      eigv(:)=vasp_eigv(:)
      kpnts(:,:)=vasp_kpnts(:,:)
      xi=vasp_xi

      ! Set up UMat2D
      tUMat2D=.TRUE.
      Allocate(UMat2D(nStates,nStates),STAT=ierr)
      Call MemAlloc(ierr,UMat2D,HElementSize*NSTATES*NSTATES,'UMAT2D')
!      call LogMemAlloc('UMat2D',nStates**2,HElementSize*8,thisroutine,tagUMat2D,ierr)
      do i=1,nStates
         do j=1,nStates
            UMat2D(i,j)=HElement(vasp_umat2d(i,j))
         end do
      end do
      
      return
   end subroutine NECIReceiveVASPData

end module vasp_interface
