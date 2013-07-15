#include  "macros.h"

! This module copies the idea of the linear scaling algorithm to use a hash function
! and index array to allow a search of the core space to be performed quickly.

module semi_stoch_hash

    use constants
    use FciMCData, only: CoreHashIndex, determ_space_size, ll_node

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
        integer :: i, ierr, DetHash, counter, total
        type(ll_node), pointer :: temp_node
        character(len=*), parameter :: t_r = "InitialiseCoreHashTable"

        allocate(CoreHashIndex(determ_space_size), stat=ierr)

        do i = 1, determ_space_size
            CoreHashIndex(i)%ind = 0
        end do

        do i = 1, determ_space_size
            call decode_bit_det(nI, core_space(:,i))
            DetHash = FindCoreHash(nI)
            temp_node => CoreHashIndex(DetHash)
            ! If the first element in the list has not been used.
            if (temp_node%ind == 0) then
                temp_node%ind = i
            else
                do while (associated(temp_node%next))
                    temp_node => temp_node%next
                end do
                allocate(temp_node%next)
                temp_node%next%ind = i
            end if
        end do

        nullify(temp_node)

    end subroutine InitialiseCoreHashTable

end module semi_stoch_hash
