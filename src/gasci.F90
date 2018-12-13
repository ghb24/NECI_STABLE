module gasci
  use SystemData, only: tNConservingGAS, tSpinConservingGAS, nBasis
  use constants
  use util_mod, only: get_free_unit
  use bit_rep_data, only: NIfTot, NIfD
  implicit none

  integer :: nGAS
  integer(n_int), allocatable :: gasOrbs(:,:)
  
  contains

    subroutine loadGAS()
      integer :: gas_unit, iOrb, nOrbs
      integer :: iGAS(1000)

      gas_unit = get_free_unit()
      open(gas_unit, file="GASOrbs",status='old')
      nOrbs = nBasis/2
      read(gas_unit,*) iGAS(1:nOrbs)
      nGAS = maxval(iGAS)
      allocate(gasOrbs(0:NIfD,nGAS))
      do iOrb = 1, nOrbs
         ! now write the orbitals read in to the current GAS
         ! set both the alpha- and the beta-spin orbital
         call setOrb(gasOrbs(:,iGAS(iOrb)),2*iOrb)
         call setOrb(gasOrbs(:,iGAS(iOrb)),2*iOrb-1)
      end do

      do iOrb = 1, nGAS
         print *, "Number of orbs in GAS", iOrb, "is", sum(popCnt(gasOrbs(:,iOrb)))
      end do

      contains 

        subroutine splitLine(line,vals,n)
          implicit none
          character(*), intent(in) :: line
          integer, intent(out) :: vals(:)
          integer, intent(out) :: n

          integer :: status, buffer(1000)
          
          n = 1
          do 
             ! use a buffer to catch exceptions
             read(line,'(A)',iostat=status) buffer(1:n)
             if(status.ne.0) exit
             vals(1:n) = buffer(1:n)
             n=n+1
          end do
          ! we overcounted by 1 because we started at 1
          n = n - 1
        end subroutine splitLine

       subroutine setOrb(ilut,orb)
          implicit none
          integer(n_int), intent(inout) :: ilut(0:NIfD)
          integer, intent(in) :: orb

          integer :: pos
          ! get the index of the integer containing this orbital
          pos = (orb-1) / bits_n_int
          ilut(pos) = ibset(ilut(pos),mod(orb-1,bits_n_int))
        end subroutine setOrb
    end subroutine loadGAS

!------------------------------------------------------------------------------------------!

    subroutine clearGAS()
      implicit none
      
      deallocate(gasOrbs)
    end subroutine clearGAS

!------------------------------------------------------------------------------------------!
    
    function isValidExcit(ilutI,ilutJ) result(valid)
      ! check if the excitation from ilutI to ilutJ is valid within the GAS
      implicit none
      integer(n_int), intent(in) :: ilutI(0:NIfTot), ilutJ(0:NIfTot)
      integer(n_int) :: gasI(0:NIfD), gasJ(0:NIfD)
      logical :: valid

      integer :: i, x(0:NIfD)
      ! integers with all even bits set
#ifdef __INT64
      integer(n_int), parameter :: evenBits = 6148914691236517205_n_int
#else
      integer(n_int), parameter :: evenBits = 1431655765_n_int
#endif

      valid = .true.
      ! safety check: do the gasOrbs exist
      if(.not. allocated(gasOrbs)) return
      do i = 1, nGAS
         gasI = gasComponent(ilutI,i)
         gasJ = gasComponent(ilutJ,i)
         x = popCnt(gasI)
         x = popCnt(gasJ)
         ! check if ilutI and ilutJ have the same number of electrons in GAS-i
         if(sum(popCnt(gasI)) .ne. sum(popCnt(gasJ))) valid = .false.
         if(tSpinConservingGAS) then
            ! check if the number of even bits set (=number of beta electrons) is the
            ! same in both GAS
            if(sum(popCnt(iand(gasI,evenBits))) .ne. sum(popCnt(iand(gasJ,evenBits)))) &
                 valid = .false.
         endif
      end do

      contains 
        
        function gasComponent(ilut, i) result(gasIlut)
          integer(n_int), intent(in) :: ilut(0:NIfTot)
          integer, intent(in) :: i
          
          integer(n_int) :: gasIlut(0:NIfD)
          
          gasIlut = iand(ilut(0:NIfD),gasOrbs(0:NIfD,i))
          
        end function gasComponent
        
    end function isValidExcit
end module gasci
