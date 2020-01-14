#include "macros.h"
module load_balance_calcnodes

    ! This module breaks circular dependencies and stores data.
    ! Otherwise it is just part of load_balance

    use FciMCData, only: HFDet, hash_iter, hash_shift
    use SystemData, only: nBasis, nel, tCSF
    use csf_data, only: csf_orbital_mask
    use MemoryManager, only: TagIntType
    use CalcData, only: tUniqueHFNode
    use bit_reps, only: NIfTot
    use Parallel_neci
    use constants

    implicit none

    integer, allocatable :: RandomOrbIndex(:), LoadBalanceMapping(:)
    integer(TagIntType) :: lb_tag
    integer :: balance_blocks
    logical :: tLoadBalanceBlocks

contains

    pure function DetermineDetNode (nel_loc, nI, iIterOffset) result(node)

        ! Depending on the Hash, determine which node determinant nI
        ! belongs to in the DirectAnnihilation scheme. NB FCIMC has each
        ! processor as a separate logical node.

        ! In:  nI   - Integer ordered list for the determinant
        ! In:  iIterOffset - Offset this iteration by this amount
        ! Out: proc - The (0-based) processor index.

        ! --> This function takes the calculated block (from get_det_block)
        !     and converts it via a simple lookup into the required node


        integer, intent(in) :: nel_loc
        integer, intent(in) :: nI(nel_loc)
        integer, intent(in) :: iIterOffset
        integer :: node

        integer :: block

        block = get_det_block(nel_loc, nI, iIterOffset)

        ! Look up the relevant node in the block-mapping.
        node = LoadBalanceMapping(block)

    end function


    pure function get_det_block(nel_loc, nI, iIterOffset) result(block)

        ! Depending on the Hash, determine which node determinant nI
        ! belongs to in the DirectAnnihilation scheme. NB FCIMC has each
        ! processor as a separate logical node.

        ! In:  nI   - Integer ordered list for the determinant
        ! In:  iIterOffset - Offset this iteration by this amount
        ! Out: proc - The (0-based) processor index.

        ! --> This function calculates the hash, and takes it mod the number
        !     of blocks in use

        integer, intent(in) :: nel_loc
        integer, intent(in) :: nI(nel_loc)
        integer, intent(in) :: iIterOffset
        integer :: block

        integer :: i
        integer(int64) :: acc
        integer(int64) :: offset
        integer(int64), parameter :: large_prime = 1099511628211_int64

        ! If we are assigning the HF to a unique processor, then put it
        ! on the last available processor.
        if (size(nI) == size(HFDet)) then
            if (tUniqueHFNode .and. all(nI == HFDet)) then
                block = nNodes
                return
            end if
        end if

        ! sum(nI) ensures that a random number is generated for each different
        ! nI, which is then added to the iteration, and the result shifted.
        ! Consequently, the less significant bits (but not the least, as these
        ! have been shifted away) for each nI will change on average once per
        ! 2^hash_shift iterations, but the change spread throughout the
        ! different iters.
        ! Generate a hash to work out an offset.  Probably very inefficient.
        if (hash_iter>0) then
           acc = 0
           do i = 1, nel_loc
               acc = (large_prime * acc) + &
                       (RandomOrbIndex(mod(iand(nI(i), csf_orbital_mask)-1,nBasis)+1) * i)
           enddo
           offset=ishft(abs(ishft(acc+hash_iter+iIterOffset, -hash_shift) ),-4)
        else
           offset=0
        endif
        acc = 0
        if(tCSF) then
            do i = 1, nel_loc
                acc = (large_prime * acc) + &
                        (RandomOrbIndex(mod(iand(nI(i), csf_orbital_mask)+offset-1,int(nBasis,int64))+1) * i)
            enddo
        else
            do i = 1, nel_loc
                acc = (large_prime * acc) + &
                     (RandomOrbIndex(mod(nI(i)+offset-1,int(nBasis,int64))+1) * i)
            enddo
        endif

        ! If the last available processor is being used for the HF det, then
        ! we can only use the the remaining (nNodes-1) processors to
        ! distribute the rest ofthe particles
        if (tUniqueHFNode) then
            block = int(abs(mod(acc, int(balance_blocks-1, int64))),sizeof_int)
        else
            block = int(abs(mod(acc, int(balance_blocks, int64))),sizeof_int)
        end if
        block = block + 1

    end function


end module
