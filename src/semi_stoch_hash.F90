#include  "macros.h"

! This module copies the idea of the linear scaling algorithm to use a hash function
! and index array to allow a search of the core space to be performed quickly.

module semi_stoch_hash

    use constants
    use FciMCData, only: CoreHashIndex, CoreHashIndexArr1, CoreHashIndexArr2, &
                         nCoreClashMax, determ_space_size

    implicit none

contains

    pure function FindCoreHash(nJ) result(hashInd)

        use csf_data, only: csf_orbital_mask
        use FciMCData, only: RandomHash
        use SystemData, only: tCSF, nel, nBasis

        integer, intent(in) :: nJ(nel)
        integer :: hashInd
        integer :: i
        integer(int64) :: hash

        hash = 0

        if(tCSF) then
            do i = 1, nel
                hash = (1099511628211_int64 * hash) + &
                        int(RandomHash(mod(iand(nJ(i), csf_orbital_mask)-1,nBasis)+1) * i,int64)
            enddo
        else
            do i = 1, nel
                hash = (1099511628211_int64 * hash) + &
                        int(RandomHash(nJ(i))*i,int64)
            enddo
        endif

        hashInd = int(abs(mod(hash, int(determ_space_size, int64))),sizeof_int)+1

    end function FindCoreHash
    
    subroutine InitialiseCoreHashTable()

        use bit_reps, only: decode_bit_det
        use FciMCData, only: core_space
        use SystemData, only: nel

        integer :: nI(nel)
        integer :: i, ierr, Slot, DetHash
        character(len=*), parameter :: t_r = "InitialiseCoreHashTable"

        nCoreClashMax = 2

        allocate(CoreHashIndexArr1(0:nCoreClashMax, determ_space_size), stat=ierr)
        if(ierr.ne.0) call stop_all(t_r, "Error in allocation")
        CoreHashIndex => CoreHashIndexArr1
        CoreHashIndex(:,:) = 0
        CoreHashIndex(0,:) = 1

        do i = 1, determ_space_size
            call decode_bit_det(nI, core_space(:,i))
            DetHash = FindCoreHash(nI)
            Slot = CoreHashIndex(0,DetHash)
            CoreHashIndex(Slot,DetHash) = i
            CoreHashIndex(0,DetHash) = CoreHashIndex(0,DetHash) + 1
            if(CoreHashIndex(0,DetHash).gt.nCoreClashMax) then
                call EnlargeCoreHashTable()
            endif
        end do

    end subroutine InitialiseCoreHashTable

    subroutine EnlargeCoreHashTable()

        integer :: ierr, i, FinalVal
        character(len=*), parameter :: t_r = "EnlargeCoreHashTable"

        nCoreClashMax = nCoreClashMax + 1

        write(6,"(a66)") "Enlarging core hash table since it is not optimally load balanced."
        write(6,"(a26,I5,a19)") "Allowing for a maximum of ",nCoreClashMax," core hash clashes."

        if(allocated(CoreHashIndexArr1)) then
            allocate(CoreHashIndexArr2(0:nCoreClashMax, determ_space_size), stat=ierr)
            if(ierr.ne.0) then
                call stop_all(t_r,"Error reallocating CoreHashIndexArr2")
            endif

            CoreHashIndexArr2(:,:)=0
            do i = 1, determ_space_size
                FinalVal = CoreHashIndex(0,i)-1
                CoreHashIndexArr2(0:FinalVal,i) = CoreHashIndex(0:FinalVal,i)
            enddo

            deallocate(CoreHashIndexArr1)
            nullify(CoreHashIndex)
            CoreHashIndex => CoreHashIndexArr2
        else
            if(.not.allocated(CoreHashIndexArr2)) then
                call stop_all(t_r,"CoreHashIndexArr2 should be allocated")
            endif
            if(allocated(CoreHashIndexArr1)) then
                call stop_all(t_r,"CoreHashIndexArr1 should not be allocated")
            endif

            allocate(CoreHashIndexArr1(0:nCoreClashMax, determ_space_size), stat=ierr)
            if(ierr.ne.0) then
                call stop_all(t_r,"Error reallocating CoreHashIndexArr2")
            endif

            CoreHashIndexArr1(:,:) = 0
            do i = 1, determ_space_size
                FinalVal = CoreHashIndex(0,i)-1
                CoreHashIndexArr1(0:FinalVal,i) = CoreHashIndex(0:FinalVal,i)
            enddo

            deallocate(CoreHashIndexArr2)
            nullify(CoreHashIndex)
            CoreHashIndex => CoreHashIndexArr1
        endif

    end subroutine EnlargeCoreHashTable

end module semi_stoch_hash
