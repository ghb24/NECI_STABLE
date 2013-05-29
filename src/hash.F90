! Copyright (c) 2013, Ali Alavi
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
#include "macros.h"

module hash

    use FciMCData , only : hash_iter, hash_shift, RandomHash
    use Parallel_neci , only : nNodes
    use constants , only : int64,sizeof_int
    use csf_data, only: csf_orbital_mask
    use Systemdata, only: nel, tCSF, nBasis

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

!sum(nI) ensures that a random number is generated for each different nI, which is then added to the iteration,
!  and the result shifted.  Consequently, the less significant bits (but not the least, as these have been shifted away)
!  for each nI will change on average once per 2^hash_shift iterations, but the change spread throughout the different iters.
        if (hash_iter>0) then  !generate a hash to work out an offset.  Probably very inefficient
           acc = 0
           do i = 1, nel
               acc = (1099511628211_int64 * acc) + &
                       (RandomHash(mod(iand(nI(i), csf_orbital_mask)-1,nBasis)+1) * i)
           enddo
           offset=ishft(abs(ishft(acc+hash_iter+iIterOffset, -hash_shift) ),-4)
        else
           offset=0
        endif
        acc = 0
        if(tCSF) then
            do i = 1, nel
                acc = (1099511628211_int64 * acc) + &
                        (RandomHash(mod(iand(nI(i), csf_orbital_mask)+offset-1,int(nBasis,int64))+1) * i)
            enddo
        else
            do i = 1, nel
                acc = (1099511628211_int64 * acc) + &
                        (RandomHash(mod(nI(i)+offset-1,int(nBasis,int64))+1) * i)
            enddo
        endif
        node = int(abs(mod(acc, int(nNodes, int64))),sizeof_int)

    end function

      
end module hash


