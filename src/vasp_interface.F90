module vasp_interface

use constants, only: dp
implicit none
public
save

integer, allocatable :: KPntInd(:)
real(dp), allocatable :: eigv(:),kpnts(:,:)
real(dp) :: xi

integer :: nKP,nStates,KPMsh(3)

contains

   subroutine NECIReceiveVASPData(vasp_nbands,vasp_nkpts,vasp_kpnts,vasp_kpmsh,vasp_xi,vasp_kpntind,vasp_eigv,vasp_umat2D,vasp_nEl)
      use constants, only: dp
      use SystemData, only: nEl
      use UMatCache, only: SetupUMatCache,UMat2D,tagUMat2D,tUMat2D
      use global_utilities
      use HElem
      implicit none
      integer :: vasp_nbands,vasp_nkpts,vasp_kpntind(vasp_nbands*vasp_nkpts),vasp_kpmsh(3)
      real(dp) :: vasp_kpnts(3,vasp_nkpts),vasp_xi,vasp_eigv(vasp_nbands*vasp_nkpts)
      integer :: vasp_nEl
#ifdef __CMPLX
      complex(dp) :: vasp_umat2d(vasp_nbands*vasp_nkpts,vasp_nbands*vasp_nkpts)
#else
      real(dp) :: vasp_umat2d(vasp_nbands*vasp_nkpts,vasp_nbands*vasp_nkpts)
#endif
      integer :: i,j,ierr
      character(*), parameter :: thisroutine='NECIReceiveVASPData'

      nKP=vasp_nkpts
      nStates=vasp_nbands*vasp_nkpts

      write (6,*) 'nkp',nKP
      write (6,*) 'kp'
      do i=1,nKP
         write (6,*) vasp_kpnts(:,i)
      end do
      nEl=vasp_nEl*nKP

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
      call LogMemAlloc('UMat2D',nStates**2,HElement_t_size*8,thisroutine,tagUMat2D,ierr)
      do i=1,nStates
         do j=1,nStates
            UMat2D(i,j)=(vasp_umat2d(i,j))
         end do
      end do

      return
   end subroutine NECIReceiveVASPData

end module vasp_interface
