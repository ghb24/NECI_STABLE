#include "macros.h"
module kMatProjE
  use SystemData, only: nBasis, tStoreSpinOrbs, tReltvy, G1, nel
  use UMatCache, only: UMatInd
  use constants
  use util_mod, only: get_free_unit, open_new_file
  use UMatCache, only: gtID
  use FciMCData, only: all_sum_proje_denominator
  use tc_three_body_data, only: kMatLin, kMatQuad, kMatLin_win, kMatQuad_win
  use shared_memory_mpi, only: shared_allocate_mpi, shared_deallocate_mpi
  use ParallelHelper, only: iProcIndex_intra
  use Parallel_neci

  implicit none

  ! their contribution to the reference energy (stored for each double excit)
  real(dp), allocatable :: kMatProjEContrib(:)
  ! size of the arrays kMat and kMatProjEContrib
  integer(int64) :: kMatSize
  logical :: tSizeSet = .false.
 
  private :: kMatProjEContrib, kMatSize, tSizeSet

  contains

!------------------------------------------------------------------------------!    

    subroutine addProjEContrib(nI,nJ,sgn)
      implicit none
      integer, intent(in) :: nI(nel),nJ(nel)
      real(dp), intent(in) :: sgn(lenof_sign)

      integer :: ex(2,2)
      logical :: tPar
      integer :: id(2,2), ind

      ! kMat reference energy contributions are coming from double excitations
      ex(1,1) = 2

      ! get the excitation from nI to nJ
      call GetExcitation(nI,nJ,nel,ex,tPar)
      id = gtID(ex)
      ! the k-matrix only couples same-spin orbitals
      if ( tReltvy.or.((G1(ex(1,1))%Ms == G1(ex(2,1))%Ms) .and. &
           (G1(ex(1,2))%Ms == G1(ex(2,2))%Ms)) ) then
         ind = UMatInd(id(1,1),id(1,2),id(2,1),id(2,2))
         ! add the average of the two spawn matrix elements (as this is what is used 
         ! in the energy calculation + dynamics)
         kmatProjEContrib(ind) =&
              kMatProjEContrib(ind) + sum(sgn)/inum_runs * &
              0.5*(kMat(ind) + kMat(UMatInd(id(1,2),id(1,1),id(2,2),id(2,1))))
      endif
      ! check the spin, do we have an exchange contribution?
      if ( tReltvy.or.((G1(ex(1,1))%Ms == G1(ex(2,2))%Ms) .and. &
           (G1(ex(1,2))%Ms == G1(Ex(2,1))%Ms)) ) then
         ind = UMatInd(id(1,1),id(1,2),id(2,2),id(2,1))
         kmatProjEContrib(ind) =&
              kMatProjEContrib(ind) - sum(sgn)/inum_runs * &
              0.5*(kMat(ind) + kMat(UMatInd(id(1,2),id(1,1),id(2,1),id(2,2))))         
      endif
    end subroutine addProjEContrib

!------------------------------------------------------------------------------!    

! print the contributions of the k-matrix to the reference energy to a file
    subroutine printProjEContrib()
      implicit none
      character(*), parameter :: refFilename = "kMatProjE"
      integer :: iunit, nBI
      integer :: i,j,k,l
      real(dp) :: tmp
      integer(MPIArg) :: ierror
      
! TODO Externalize the spinorb check
      if(tStoreSpinOrbs) then
         nBI = nBasis
      else
         nBI = nBasis/2
      endif

! I/O only done by root
      if(iProcIndex.eq.root) then
         ! do an in-place reduction to save memory
         call MPI_Reduce(MPI_IN_PLACE,kMatProjEContrib,kMatSize,MPI_DOUBLE_PRECISION,&
              root,MPI_SUM,MPI_COMM_WORLD,ierror)
         ! open the file
         iunit = get_free_unit()
         call open_new_file(iunit, refFilename)

         ! for each possible double excit, print the contribution
         ! TODO: In principle, we only need to loop over i,j occ and k,l virtual
         do i = 1, nBI
            do j = i, nBI
               do k = 1, nBI
                  do l = k, nBI
                     tmp = kMatProjEContrib(UMatInd(i,j,k,l))
                     if(abs(tmp)>eps) then
                        ! stored is the full accumulated contribution, but we are
                        ! only interested in the average
                        tmp = tmp / sum(all_sum_proje_denominator)
                        write(iunit,*) i,j,k,l,tmp
                     endif
                  end do
               end do
            end do
         end do

         tmp = sum(kMatProjEContrib)/sum(all_sum_proje_denominator)
         write(iout,*) "Total transcorrelated 2-Body correlation energy", tmp
      else
         call MPI_Reduce(kMatProjEContrib,kMatProjEContrib,kMatSize,MPI_DOUBLE_PRECISION,&
              root,MPI_SUM,MPI_COMM_WORLD,ierror)
      endif
    end subroutine printProjEContrib

!------------------------------------------------------------------------------!    

    subroutine readKMat()
      implicit none

      call determineKMatSize()
      call readKMatComponent(kMatLin_win, kMatLin, "KDUMPLIN")
      call readKMatComponent(kMatQuad_win, KMatQuad, "KDUMPQUAD")

    end subroutine readKMat

!------------------------------------------------------------------------------!    

! read in the transcorrelated part of the two-body integrals
    subroutine readKMatComponent(kMatArr_win, kMatArr, filename)
      implicit none
      real(dp), pointer :: kMatArr(:)
      integer(MPIArg), intent(in) :: kMatArr_win
      character(*) :: filename
      integer :: iunit, ierr
      integer :: i,j,k,l
      real(dp) :: matel
      character(*), parameter :: t_r = "readKMatComponent"
      
      ! allocate the containers
      call setupKMat(kMatArr_win, kMatArr)

      ! only have root read per node
      if(iProcIndex_intra.eq.0) then
         ! open the file
         iunit = get_free_unit()
         open(iunit,file = filename, status = 'old')
         ! read the integrals
         do 
            read(iunit,*,iostat = ierr) matel, i, k, j, l
            if(ierr < 0) then
               exit
            else if(ierr > 0) then
               call stop_all(t_r, "Error reading KDUMP file")
            else
               ! store the matrix element
               kMatArr(UMatInd(i,j,k,l)) = matel
            end if
         end do
      endif
    end subroutine readKMatComponent

!------------------------------------------------------------------------------!

    subroutine determineKMatSize()
      implicit none
      integer :: nBI

      if(tStoreSpinOrbs) then
         nBI = nBasis
      else
         nBI = nBasis / 2
      end if

! kmat has the same access pattern as umat, so use UMatInd as indexing function
! the index of the largest element
      kMatSize = UMatInd(nBI,nBI,nBI,nBI)
      tSizeSet = .true.
    end subroutine determineKMatSize

!------------------------------------------------------------------------------!    

! allocate the container for the transcorrelated two-body terms
! and the container for the reference energy contributions
    subroutine setupKMat(kMatArr_win, kMatArr)
      implicit none
      real(dp), pointer :: kMatArr(:)
      integer(MPIArg), intent(in) :: kMatArr_win
      integer :: ierr
      character(*), parameter :: t_r = "setupKMat"

      call shared_allocate_mpi(kmatArr_win, kMatArr, (/kMatSize/)) 
      kMatArr = 0.0_dp
      allocate(kMatProjEContrib(kMatSize), stat = ierr)
      kMatProjEContrib = 0.0_dp
      if(ierr.ne.0) call stop_all(t_r,"Allocation failed")

    end subroutine setupKMat

!------------------------------------------------------------------------------!    

    subroutine freeKMat()
      implicit none
      
      if(associated(kMatLin)) call shared_deallocate_mpi(kMatLin_win, kMatLin)
      if(associated(kMatQuad)) call shared_deallocate_mpi(kMatQuad_win, kMatQuad)
      if(allocated(kMatProjEContrib)) deallocate(kMatProjEContrib)
    end subroutine freeKMat

!------------------------------------------------------------------------------!    

    function kMat(index) result(val)
      implicit none
      integer(int64), intent(in) :: index
      real(dp) :: val
      
      val = kMatLin(index) + kMatQuad(index)
    end function kMat

!------------------------------------------------------------------------------!    
    
end module kMatProjE
