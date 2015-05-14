#include "macros.h"
module load_balance

    use FciMCData, only: HFDet, hash_iter, hash_shift
    use csf_data, only: csf_orbital_mask
    use CalcData, only: tUniqueHFNode
    use Parallel_neci, only: nNodes
    use SystemData, only: tCSF, nBasis
    use constants
    use util_mod

    implicit none

    integer, allocatable :: RandomOrbIndex(:), LoadBalanceMapping(:)
    integer(TagIntType) :: lb_tag
    integer :: balance_blocks
    logical :: tLoadBalanceBlocks

    ! TODO:
    ! - Initialise mapping
    ! - Modify DetermineDetNode to use the mapping
    ! - Add mechanism to shuffle walkers around
    ! - Add redistribution of walkers
    ! - Integrate with POPSFILES. Need to output the mapping for restarts.

contains

    subroutine init_load_balance()

        ! Initialise the load balancing.
        !
        ! n.b. The initialisation of RandomOrbIndex remains in SetupParameters
        !      to preserve sequencing, which maintains testcode results.

        integer :: oversample_factor, ierr, i
        character(*), parameter :: t_r = 'init_load_balance'

        !
        ! Initialise the mapping of balancing blocks to nodes. By default
        ! this is just a uniform mapping.
        ASSERT(.not. (tLoadBalanceBlocks .and. tUniqueHFNode))
        if (tLoadBalanceBlocks) then
            oversample_factor = 100
        else
            oversample_factor = 1
        end if

        balance_blocks = nNodes * oversample_factor
        allocate(LoadBalanceMapping(balance_blocks), stat=ierr)
        log_alloc(LoadBalanceMapping, lb_tag, ierr)

        ! Generate a uniform mapping(by default)
        do i = 1, balance_blocks
            LoadBalanceMapping(i) = int((i - 1) / oversample_factor)
        end do

    end subroutine

    subroutine clean_load_balance()

        character(*), parameter :: this_routine = 'clean_load_balance'

        if (allocated(LoadBalanceMapping)) then
            deallocate(LoadBalanceMapping)
            log_dealloc(lb_tag)
        end if

    end subroutine

    pure function DetermineDetNode (nel_loc, nI, iIterOffset) result(node)

        ! Depending on the Hash, determine which node determinant nI
        ! belongs to in the DirectAnnihilation scheme. NB FCIMC has each
        ! processor as a separate logical node.

        ! In:  nI   - Integer ordered list for the determinant
        ! In:  iIterOffset - Offset this iteration by this amount
        ! Out: proc - The (0-based) processor index.

        integer, intent(in) :: nel_loc
        integer, intent(in) :: nI(nel_loc)
        integer, intent(in) :: iIterOffset
        integer :: node
        
        integer :: i, block
        integer(int64) :: acc
        integer(int64) :: offset
        integer(int64), parameter :: large_prime = 1099511628211_int64

        ! If we are assigning the HF to a unique processor, then put it 
        ! on the last available processor.
        if (size(nI) == size(HFDet)) then
            if (tUniqueHFNode .and. all(nI == HFDet)) then
                node = nNodes-1
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

        ! Look up the relevant node in the block-mapping.
        node = LoadBalanceMapping(block)

    end function DetermineDetNode


end module
