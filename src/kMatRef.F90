#include "macros.h"
module kMatRef
  use SystemData, only: nBasis, tStoreSpinOrbs
  use UMatCache, only: UMatInd
  use constants
  use util_mod, only: get_free_unit

  implicit none

  real(dp), allocatable :: kMat(:)
  real(dp), allocatable :: kMatRefContrib(:)

  private :: kMat, kMatRefContrib

  contains

!------------------------------------------------------------------------------!    

    subroutine addRefContrib(i,j,k,l,sgn)
      implicit none
      integer, intent(in) :: i,j,k,l
      real(dp), intent(in) :: sgn(lenof_sign)

      ! add the average of the two spawn matrix elements (as this is what is used 
      ! in the energy calculation + dynamics)
      kmatRefContrib(UMatInd(i,j,k,l)) = kMatRefContrib(UMatInd(i,j,k,l)) + sum(sgn)/inum_runs * &
           0.5*(kMat(UMatInd(i,j,k,l)) + kMat(UMatInd(j,i,l,k)))
    end subroutine addRefContrib

!------------------------------------------------------------------------------!    

! read in the transcorrelated part of the two-body integrals
    subroutine readKMat()
      implicit none
      integer :: iunit, ierr
      integer :: i,j,k,l
      real(dp) :: matel
      character(*), parameter :: kMatFilename = "KDUMP"
      character(*), parameter :: t_r = "readKMat"
      
      ! allocate the containers
      call setupKMat()

      ! open the file
      iunit = get_free_unit()
      open(iunit,file = kMatFilename, status = 'old')
      ! read the integrals
      do 
         read(iunit,*,iostat = ierr) matel, i, k, j, l
         if(ierr < 0) then
            exit
         else if(ierr > 0) then
            call stop_all(t_r, "Error reading KDUMP file")
         else
            ! store the matrix element
            kMat(UMatInd(i,j,k,l)) = matel
         end if
      end do
    end subroutine readKMat

!------------------------------------------------------------------------------!    

! allocate the container for the transcorrelated two-body terms
! and the container for the reference energy contributions
    subroutine setupKMat()
      implicit none
      integer :: size, nBI, ierr
      character(*), parameter :: t_r = "setupKMat"

      if(tStoreSpinOrbs) then
         nBI = nBasis
      else
         nBI = nBasis / 2
      end if

! kmat has the same access pattern as umat, so use UMatInd as indexing function
! the index of the largest element
      size = UMatInd(nBI,nBI,nBI,nBI)

      allocate(kMat(size), stat = ierr)
      kMat = 0.0_dp
      if(ierr.ne.0) call stop_all(t_r,"Allocation failed")
      allocate(kMatRefContrib(size), stat = ierr)
      kMatRefContrib = 0.0_dp
      if(ierr.ne.0) call stop_all(t_r,"Allocation failed")

    end subroutine setupKMat

!------------------------------------------------------------------------------!    

    subroutine clearKMat()
      implicit none
      
      if(allocated(kMat)) deallocate(kMat)
      if(allocated(kMatRefContrib)) deallocate(kMatRefContrib)
    end subroutine clearKMat

!------------------------------------------------------------------------------!    
    
end module kMatRef
