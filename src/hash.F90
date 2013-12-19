#include "macros.h"

module hash

    use bit_reps, only: set_flag
    use FciMCData , only : hash_iter, hash_shift, RandomHash, RandomHash2, HFDet
    use Parallel_neci , only : nNodes
    use constants , only : int64,sizeof_int
    use csf_data, only: csf_orbital_mask
    use Systemdata, only: nel, tCSF, nBasis
    use CalcData, only: tUniqueHFNode

    implicit none

    contains

    pure function DetermineDetNode (nI, iIterOffset) result(node)

        ! Depending on the Hash, determine which node determinant nI
        ! belongs to in the DirectAnnihilation scheme. NB FCIMC has each processor as a separate logical node.
        !
        ! In:  nI   - Integer ordered list for the determinant
        ! In:  iIterOffset - Offset this iteration by this amount
        ! Out: proc - The (0-based) processor index.

        implicit none
        integer, intent(in) :: nI(nel)
        integer, intent(in) :: iIterOffset
        integer :: node
        
        integer :: i
        integer(int64) :: acc
        integer(int64) :: offset
        integer(int64), parameter :: large_prime = 1099511628211_int64

        ! If we are assigning the HF to a unique processor, then put it 
        ! on the last available processor.
        if (tUniqueHFNode .and. all(nI == HFDet)) then
            node = nNodes-1
            return
        end if

!sum(nI) ensures that a random number is generated for each different nI, which is then added to the iteration,
!  and the result shifted.  Consequently, the less significant bits (but not the least, as these have been shifted away)
!  for each nI will change on average once per 2^hash_shift iterations, but the change spread throughout the different iters.
        if (hash_iter>0) then  !generate a hash to work out an offset.  Probably very inefficient
           acc = 0
           do i = 1, nel
               acc = (large_prime * acc) + &
                       (RandomHash(mod(iand(nI(i), csf_orbital_mask)-1,nBasis)+1) * i)
           enddo
           offset=ishft(abs(ishft(acc+hash_iter+iIterOffset, -hash_shift) ),-4)
        else
           offset=0
        endif
        acc = 0
        if(tCSF) then
            do i = 1, nel
                acc = (large_prime * acc) + &
                        (RandomHash(mod(iand(nI(i), csf_orbital_mask)+offset-1,int(nBasis,int64))+1) * i)
            enddo
        else
            do i = 1, nel
                acc = (large_prime * acc) + &
                        (RandomHash(mod(nI(i)+offset-1,int(nBasis,int64))+1) * i)
            enddo
        endif

        ! If the last available processor is being used for the HF det, then
        ! we can only use the the remaining (nNodes-1) processors to
        ! distribute the rest ofthe particles
        if (tUniqueHFNode) then
            node = int(abs(mod(acc, int(nNodes-1, int64))),sizeof_int)
        else
            node = int(abs(mod(acc, int(nNodes, int64))),sizeof_int)
        end if

    end function DetermineDetNode

    ! Routin to find the correct position in the hash table
    pure function FindWalkerHash(nJ, HashIndexLength) result(hashInd)
        implicit none
        integer, intent(in) :: nJ(nel)
        integer, intent(in) :: HashIndexLength
        integer :: hashInd
        integer :: i
        integer(int64) :: hash
        hash = 0
!        write(6,*) "nJ: ",nJ(:)
        if(tCSF) then
            do i = 1, nel
                hash = (1099511628211_int64 * hash) + &
                        int(RandomHash2(mod(iand(nJ(i), csf_orbital_mask)-1,nBasis)+1) * i,int64)
            enddo
        else
            do i = 1, nel
                hash = (1099511628211_int64 * hash) + &
                        int(RandomHash2(nJ(i))*i,int64)

!                        (RandomHash(mod(nI(i)+offset-1,int(nBasis,int64))+1) * i)
            enddo
        endif
        hashInd = int(abs(mod(hash, int(HashIndexLength, int64))),sizeof_int)+1

    end function FindWalkerHash

    FUNCTION CreateHash(DetCurr)
        INTEGER :: DetCurr(NEl),i
        INTEGER(KIND=int64) :: CreateHash

        CreateHash=0
        do i=1,NEl
!            CreateHash=13*CreateHash+i*DetCurr(i)
            CreateHash=(1099511628211_int64*CreateHash)+i*DetCurr(i)
            
!            CreateHash=mod(1099511628211*CreateHash,2**64)
!            CreateHash=XOR(CreateHash,DetCurr(i))
        enddo
!        WRITE(6,*) CreateHash
        RETURN

    END FUNCTION CreateHash

    subroutine reset_hash_table(hash_table)

        type(ll_node), pointer, intent(in) :: hash_table(:)
        type(ll_node), pointer :: curr, prev

        ! Reset the hash index array.
        do i = 1, size(hash_table)
            curr => hash_table(i)%next
            prev => hash_table(i)
            prev%ind = 0
            nullify(prev%next)
            do while (associated(curr))
                prev => curr
                curr => curr%next
                deallocate(prev)
            end do
        end do
        nullify(curr)
        nullify(prev)

    end subroutine reset_hash_table
      
end module hash


