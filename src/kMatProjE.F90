#include "macros.h"
module kMatProjE
  use SystemData, only: nBasis, tStoreSpinOrbs, tReltvy, G1, nel, nBI
  use UMatCache, only: UMatInd
  use constants
  use util_mod, only: get_free_unit, open_new_file
  use UMatCache, only: gtID
  use FciMCData, only: all_sum_proje_denominator
  use shared_memory_mpi, only: shared_allocate_mpi, shared_deallocate_mpi
  use ParallelHelper, only: iProcIndex_intra
  use Parallel_neci
  use LoggingData, only: tLogKMatProjE

  implicit none

  ! their contribution to the reference energy (stored for each double excit)
  real(dp), allocatable :: kMatProjEContrib(:)
  real(dp), parameter :: kMatLinFac = 0.25
  real(dp), parameter :: kMatSqFac = 0.375
  real(dp), parameter :: kMatParFac = 0.5
 
  private :: kMatProjEContrib

  type :: kMat_t
     private
     ! mpi shared memory window
     integer(MPIArg) :: shm_win
     ! pointer to the allocated array
     real(dp), pointer :: kMat_p(:)
     ! size of the array
     integer(int64) :: kMatSize

     ! member functions
     contains
       ! initialization routines
       procedure, public :: readKMatFromFile
       procedure, public :: setupKMat
       ! finalization routine (should be a destructor)
       procedure, public :: freeMemory
       ! exchange/direct matrix elements
       procedure, public :: directElement
       procedure, public :: exchElement

       ! getter for elements of kMat_p
       procedure, public :: elementAccess
  end type kMat_t

  type(kMat_t) :: kMatLin, kMatSq
  type(kMat_t) :: kMatAA, kMatAB

  contains


!------------------------------------------------------------------------------!  

    subroutine freeMemory(this) 
      implicit none
      class(kMat_t) :: this
      
      if(associated(this%kMat_p)) call shared_deallocate_mpi(this%shm_win, this%kMat_p)
    end subroutine freeMemory

!------------------------------------------------------------------------------!    

    function directElement(this, i,j,k,l) result(matel)
      implicit none
      class(kMat_t) :: this
      integer, intent(in) :: i,j,k,l
      real(dp) :: matel
      
      matel = this%kMat_p(UMatInd(i,j,k,l))
    end function directElement

!------------------------------------------------------------------------------!    

    function exchElement(this, i,j,k,l) result(matel)
      implicit none
      class(kMat_t) :: this
      integer, intent(in) :: i,j,k,l
      real(dp) :: matel
      
      matel = this%kMat_p(UMatInd(i,j,l,k))
    end function exchElement

!------------------------------------------------------------------------------!      

! read in the transcorrelated part of the two-body integrals
    subroutine readKMatFromFile(this, filename)
      implicit none
      class(kMat_t) :: this
      character(*) :: filename
      integer :: iunit, ierr
      integer :: i,j,k,l
      real(dp) :: matel
      character(*), parameter :: t_r = "readKMatFromFile"
      
      ! allocate the containers
      call this%setupKMat()

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
               this%kMat_p(UMatInd(i,j,k,l)) = matel
            end if
         end do
      endif
    end subroutine readKMatFromFile

!------------------------------------------------------------------------------!    

! allocate the container for the transcorrelated two-body terms
! and the container for the reference energy contributions
    subroutine setupKMat(this)
      implicit none
      class(kMat_t) :: this
      character(*), parameter :: t_r = "setupKMat"

      this%kMatSize = determineKMatSize()
      call shared_allocate_mpi(this%shm_win, this%kMat_p, (/this%kMatSize/)) 
      this%kMat_p = 0.0_dp
    end subroutine setupKMat

!------------------------------------------------------------------------------!    

    subroutine freeKMat()
      implicit none
      
      call kMatLin%freeMemory()
      call kMatSq%freeMemory()
      if(allocated(kMatProjEContrib)) deallocate(kMatProjEContrib)
    end subroutine freeKMat

!------------------------------------------------------------------------------!    
! Specific member functions for direct element access
!------------------------------------------------------------------------------!    

    function elementAccess(this, index) result(kMatel)
      implicit none
      class(kMat_t) :: this
      integer(int64), intent(in) :: index
      real(dp) :: kMatel
      
      kMatel = this%kMat_p(index)
    end function elementAccess

!------------------------------------------------------------------------------!    
! Generic non-member functions used for the projected energy estimate
!------------------------------------------------------------------------------!    

    subroutine addProjEContrib(nI,nJ,sgn)
      implicit none
      integer, intent(in) :: nI(nel),nJ(nel)
      real(dp), intent(in) :: sgn(lenof_sign)

      integer :: ex(2,2)
      logical :: tPar
      integer :: id(2,2)
      integer(int64) :: ind

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
      integer :: iunit
      integer :: i,j,k,l
      real(dp) :: tmp
      integer(MPIArg) :: ierror
      integer :: kMatSize
      
      kMatSize = determineKMatSize()

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
      character(*), parameter :: t_r = "readKMat"
      integer :: ierr, kMatSize
      call kMatLin%readKMatFromFile("KDUMPLIN")
      call kMatSq%readKMatFromFile("KDUMPSQ")

      ! now, take care of the projected energy contribution if required
      if(tLogKMatProjE) then
         kMatSize = determineKMatSize()
         allocate(kMatProjEContrib(kMatSize), stat = ierr)
         kMatProjEContrib = 0.0_dp
         if(ierr.ne.0) call stop_all(t_r,"Allocation failed")
      endif

    end subroutine readKMat

!------------------------------------------------------------------------------!    

    subroutine readSpinKMat()
      implicit none
      
      call kMatAA%readKMatFromFile("KDUMPAA")
      call kMatAB%readKMatFromFile("KDUMPAB")
    end subroutine readSpinKMat

!------------------------------------------------------------------------------!

    function determineKMatSize() result(kMatSize)
      implicit none
      integer(int64) :: kMatSize

! kmat has the same access pattern as umat, so use UMatInd as indexing function
! the index of the largest element
      kMatSize = UMatInd(int(nBI),int(nBI),int(nBI),int(nBI))
    end function determineKMatSize

!------------------------------------------------------------------------------!    

    function kMat(index) result(val)
      implicit none
      integer(int64), intent(in) :: index
      real(dp) :: val
      
      val = kMatLin%elementAccess(index) + kMatSq%elementAccess(index)
    end function kMat

!------------------------------------------------------------------------------!    

    function kMatOppSpinCorrection(i,j,k,l) result(matel)
      implicit none
      integer, intent(in) :: i,j,k,l
      real(dp) :: matel

      matel = kMatLinFac * kMatLin%directElement(i,j,k,l) &
           + kMatSqFac * kMatSq%directElement(i,j,k,l) &
           - kMatLinFac * kMatLin%exchElement(i,j,k,l) &
           - kMatSqFac * kMatSq%exchElement(i,j,k,l)
    end function kMatOppSpinCorrection

!------------------------------------------------------------------------------!    

    function kMatParSpinCorrection(i,j,k,l) result(matel)
      implicit none
      integer, intent(in) :: i,j,k,l
      real(dp) :: matel

      matel = kMatParFac * kMat(UMatInd(i,j,k,l))
    end function kMatParSpinCorrection

!------------------------------------------------------------------------------!    

    function spinKMatContrib(i,j,k,l,s1,s2) result(matel)
      implicit none
      integer, intent(in) :: i,j,k,l,s1,s2
      real(dp) :: matel

      ! kMat enters the Hamiltonian with a negative sign - add here
      if(s1.eq.s2) then
         matel = - kMatAA%directElement(i,j,k,l)
      else
         matel = - kMatAB%directElement(i,j,k,l)
      endif
    end function spinKMatContrib

!------------------------------------------------------------------------------!    


    
end module kMatProjE
