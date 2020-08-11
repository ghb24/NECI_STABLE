#include "macros.h"

module LMat_freeze
    use IntegralsData, only: nFrozen, UMat
    use UMatCache, only: numBasisIndices, UMatInd
    use constants
    use util_mod, only: operator(.div.)
    use OneEInts, only: TMat2D
    use SystemData, only: ECore
    implicit none

    private
    public :: freeze_lmat, t_freeze, map_indices

    !> Parameter for number of double excitations a single LMat entry can contribute to
    !! (counting permutations)
    integer, parameter :: num_ex = 8
    integer, parameter :: num_inds = 6
    ! Step corresponding to conjugation of an index
    integer, parameter :: step = (num_inds / 2)

    ! Direct orbital excitation = pair of indices in an index array which is apart by step
        
contains

    !------------------------------------------------------------------------------------------!

    !> Checks if an entry is zeroed due to frozen orbitals being included
    !> @param[in] indices  array of size 6 containing the indices of the entry in question in the frozen orbital numbering
    !> @return t_freeze  true if the entry is zeroed
    function t_freeze(indices)
        integer(int64), intent(in) :: indices(num_inds)
        logical :: t_freeze

        t_freeze = any(indices < 1)
    end function t_freeze

    !------------------------------------------------------------------------------------------!  

    !> Maps a set of six indices from pre-freeze to post-freeze orbital indexing
    !> @param indices  on entry: array of indices in pre-freeze indexing, on return: same array in post-freeze indexing
    subroutine map_indices(indices)
        integer(int64), intent(inout) :: indices(:)

        indices = indices - numBasisIndices(nFrozen)
    end subroutine map_indices

    !------------------------------------------------------------------------------------------!  

    !> Checks if the entry is neglected due to frozen orbitals being included and absorbs entries
    !! into the lower order matrix elements if required
    subroutine freeze_lmat(matel, indices)
        HElement_t(dp), intent(inout) :: matel
        integer(int64), intent(inout) :: indices(num_inds)

        ! Offset the orbital indexing
        call map_indices(indices)
        call add_core_en(matel, indices)
    end subroutine freeze_lmat

    !------------------------------------------------------------------------------------------!    

    !> Absorb entries with repeated frozen orbitals into the corresponding lower-order
    !! terms.
    subroutine add_core_en(matel, indices)
        HElement_t(dp), intent(inout) :: matel
        integer(int64), intent(in) :: indices(num_inds)

        integer :: counts
        integer(int64) :: idx(num_ex), ct
        real(dp) :: prefactor(num_ex)

        if(t_freeze(indices)) then
            ! Count the number of different indices appearing
            counts = count_frozen_inds(indices)
            ! How many pairs of frozen orbitals do we have
            ! 1 => double excitation
            if(counts == 1) then 
                idx = frozen_double_entry(indices, prefactor)
                ! There are up to eight potential double excitations with four given indices
                ! to which the LMat entry contributes - this is because LMat is hermitian
                ! but UMat is not
                do ct = 1, num_ex
                    ! If the index is assigned, there is a duplicate frozen orb
                    if(idx(ct) > 0) then
                        ! Absorb the matrix element into UMat (with the corresponding prefactor
                        ! +/- 1/2 according to permutation and number of possible spin terms
                        UMat(idx(ct)) = UMat(idx(ct)) + prefactor(ct)*matel
                    endif
                end do
                ! 2 => single excitation
            elseif(counts == 2) then
                idx(1:2) = frozen_single_entry(indices, prefactor(1))
                
                if(idx(1) > 0) then
                    ! Absorbt the matrix element into TMat
                    ! TMat2D is indexed with spin orbs
                    call add_to_tmat(spatToSpinAlpha(idx(1)),spatToSpinAlpha(idx(2)))
                    call add_to_tmat(spatToSpinBeta(idx(1)),spatToSpinBeta(idx(2)))
                endif
                ! 3 => diagonal element
            elseif(counts == 3) then
                call frozen_diagonal_entry(indices, prefactor(1))
                ECore = ECore + prefactor(1) * matel
            endif
            ! Zero the matrix element for further usage (i.e. will not turn up anymore)
            matel = 0.0_dp
        endif

    contains
        
        subroutine add_to_tmat(ind1, ind2)
            integer(int64), intent(in) :: ind1, ind2
            
            TMat2D(ind1,ind2) = TMat2D(ind1,ind2) + prefactor(1) * matel
        end subroutine add_to_tmat

    end subroutine add_core_en

    !------------------------------------------------------------------------------------------!    

    pure function count_frozen_inds(indices) result(level)
        integer(int64), intent(in) :: indices(num_inds)
        integer :: level
        integer :: counts

        integer :: ct
        level= 0
        ! If there is an unpaired frozen orbital, the entry is discarded, it does not
        ! contribute to any contraction
        do ct = 1, num_inds
            if(count(indices(ct) == indices) == 1 .and. indices(ct) < 1) return
        end do
        ! Else, it contributes to a contraction according to the number of frozen indices
        counts = count(indices < 1)
        ! An even number of frozen orbital indices is required
        if(modulo(counts,2) /= 0) return

        level = counts .div. 2

    end function count_frozen_inds

    !------------------------------------------------------------------------------------------!    

    !> Get the index of the UMat entry to which the LMat entry with given indices shall be
    !! added if it is frozen.
    !> @param[in] indices  array of lenght 6 with the orbital indices of the LMat entry
    !> @param[ou] t_par  flag indicating if the matrix element enters UMat with a -1
    !> @return index  index of the UMat entry to add the LMat entry to, 0 if entry is not frozen
    function frozen_double_entry(indices, prefactor) result(idx)
        integer(int64), intent(in) :: indices(num_inds)
        real(dp), intent(out) :: prefactor(num_ex)
        integer(int64) :: idx(num_ex)
        integer :: f_one, f_two
        logical :: unfrozen(num_inds), t_per
        integer(int64) :: uf_idx(4)

        idx = 0        
        unfrozen = indices > 0
        ! The indices of the UMat entries are given by the non-frozen indices
        uf_idx = pack(indices,unfrozen)
        idx = permute_umat_inds(int(uf_idx(1)),int(uf_idx(2)),int(uf_idx(3)),int(uf_idx(4)),prefactor)

        ! Check the permutation and possible factor of two due to spin
        ! First, get the two frozen orbitals
        f_one = findloc(unfrozen, .false., dim = 1, back = .false.)
        f_two = findloc(unfrozen, .false., dim = 1, back = .true.)

        ! If an odd number of steps is required to end up with a direct excitation, 
        ! a factor of -1 is added
        t_per = .false.
        ! This is the case if a) both indices are src/tgt and have an odd distance
        if(f_one > step .or. f_two <= step) then
            if(modulo(f_two - f_one,2) == 1) t_per = .true.
        ! Or there is one upper and one lower index with an even distance
        elseif(modulo(f_two - f_one,2) == 0) then
            t_per = .true.
        endif
        
        if(t_per) prefactor = -1.0_dp * prefactor
        ! There are two spin configurations in a close-shell frozen scenario for terms
        ! with a direct frozen orbital (i.e. no permutation required)        
        if(is_direct(f_one,f_two)) prefactor = 2.0_dp * prefactor

    end function frozen_double_entry

    !------------------------------------------------------------------------------------------!    

    !> Returns the UMatInd values of all possible permutations of the input indices
    function permute_umat_inds(a,b,c,d,prefactor) result(inds)
        integer, intent(in) :: a,b,c,d
        real(dp), intent(out) :: prefactor(num_ex)
        integer(int64) :: inds(num_ex)

        integer :: ct

        inds = 0
        ! These are adjoint to each other, they are the same in LMat, but not UMat
        inds(1) = UMatInd(a,b,c,d)
        inds(2) = UMatInd(c,b,a,d)
        inds(3) = UMatInd(a,d,c,b)
        inds(4) = UMatInd(c,d,a,b)
        ! They have a prefactor of +1
        prefactor(1:4) = 1.0_dp
        ! These are the exchange terms, they stem from the LMat exchange terms and enter with -1
        inds(5) = UMatInd(a,b,d,c)
        inds(6) = UMatInd(d,b,a,c)
        inds(7) = UMatInd(a,c,d,b)
        inds(8) = UMatInd(d,c,a,b)
        prefactor(5:8) = -1.0_dp

        ! All terms are only counted if they are different from a previously occuring index
        ! This approach might be less performant than competing ways of doing this check,
        ! but it is very clear and less prone to bugs. If this proves to be a performance
        ! bottleneck, it can be changed to a faster implementation later on
        do ct = 2, num_ex
            ! If the index already appeared, delete it
            if(any(inds(ct) == inds(1:ct-1))) inds(ct) = 0
        end do
    end function permute_umat_inds

    !------------------------------------------------------------------------------------------!

    function frozen_single_entry(indices, prefactor) result(orbs)
        integer(int64), intent(in) :: indices(num_inds)
        real(dp), intent(out) :: prefactor

        integer :: orbs(2)
        integer :: f_orbs(num_inds-2), ct
        integer :: uf_one, uf_two
        logical :: unfrozen(num_inds)

        orbs = 0
        unfrozen = indices > 0
        f_orbs = pack(indices, .not. unfrozen)
        ! For the single entry, there are only two unfrozen orbs, find these and return them
        orbs = pack(indices, unfrozen)

        ! Again, two things to check: The number of direct orbitals and the permutation

        ! Get the posisitons of the unfrozen orbs
        uf_one = findloc(unfrozen, .true., dim = 1, back = .false.)
        uf_two = findloc(unfrozen, .true., dim = 1, back = .true.)

        ! At this point, there is no unpaired frozen orbital, so there are two options:
        ! a) there are two different orbs
        prefactor = 1.0_dp
        if(count(f_orbs(1) == f_orbs) < size(f_orbs)) then
            ! We have to count how many potential spin configurations contribute

            ! The parity is counting the number of permutations required to end up
            ! with a completely direct expression
            ! -> The parity is odd if and only if the number of direct excitations is exactly one
            ! (then, one permutation is required to fix the two non-direct excits, else,
            ! (either zero or two are needed)
            ! -> count every direct excitation with a factor -1
            do ct = 1, step
                ! Each direct repeated orbital doubles the number of spin configs
                ! (the case of four identical orbs is excluded)                
                if(is_direct_pair(indices,ct) .and. .not. unfrozen(ct)) then
                    prefactor = -2.0_dp * prefactor
                    ! Only do this once, the second occurence coincides with a direct
                    ! unfrozen orbital excitation and is thus covered below
                    exit
                endif
            end do
            ! If the single excitation is direct, its spin does not relate to the
            ! spin of the frozen orbitals -> both configs contribute
            if(is_direct(uf_one,uf_two)) prefactor = -2.0_dp * prefactor
        else
        ! b) All frozen orbs are the same -> only one spin config is allowed, only parity of
        ! the permutation counts
            ! Here, only a direct exctiation has even parity (same rule as above)
            if(.not. is_direct(uf_one, uf_two)) prefactor = -1.0_dp * prefactor
        end if
        
    end function frozen_single_entry

    !------------------------------------------------------------------------------------------!    
    
    !> Determine the prefactor for a diagonal contribution (i.e. three pairs of frozen orbitals)
    !> @param[in] indices  orbital indices of the frozen entry - need to belong to a contribution
    !                      to a diagonal matrix element of only frozen orbitals
    !> @param[out] prefactor  prefactor with which this LMat entry enters the matrix element
    pure subroutine frozen_diagonal_entry(indices, prefactor)
        integer(int64), intent(in) :: indices(num_inds)
        real(dp), intent(out) :: prefactor

        integer :: ct, directs
        logical :: t_quad
        ! Only one info is needed here, that is the parity and number of contributing spin
        ! configs. 

        ! If there are five or more repeated indices, the LMat entry does not contribute
        ! (it would require three same-spin electrons
        do ct = 1,3
            if(count(indices(ct) == indices) > 4) then
                prefactor = 0.0_dp
                return
            end if
            if(count(indices(ct) == indices) > 3) t_quad = .true.
        end do

        ! Both relate to the number of different direct excits, so count these
        directs = 0
        do ct = 1, step
            if(is_direct_pair(indices, ct)) directs = directs + 1
        end do
        
        ! If there is a quadruple index, the spin of it is fixed (both have to be occupied), then,
        ! the prefactor is always +/- 2, depending on if there is more than one direct pair
        if(t_quad) then
            prefactor = 2.0_dp
            ! The only options with a quadruple index are
            ! a) one direct excitation => -2
            ! b) three direct exctiations => +2
            if(directs == 1) prefactor = -1.0_dp * prefactor
        else
            ! In the other case, there are three relevant cases:
            ! a) All excitations are direct => prefactor 8 = 2**3 (all spins are free)
            ! b) One excitation is direct => prefactor -4 (from the spin freedom of the direct orb + overall spin)
            ! c) No excitation is direct => prefactor 2 (overall spin)
            if(directs == 0) then
                prefactor = 2.0_dp
            elseif(directs == 1) then
                prefactor = -4.0_dp
            else ! directs == 3
                prefactor = 8.0_dp
            endif
        endif
        
    end subroutine frozen_diagonal_entry

    !------------------------------------------------------------------------------------------!

    pure function is_direct(one, two) result(t_dir)
        integer, intent(in) :: one, two
        logical :: t_dir

        t_dir = (one + step == two)
    end function is_direct

    !------------------------------------------------------------------------------------------!

    pure function is_direct_pair(indices,ct) result(t_dir)
        integer(int64), intent(in) :: indices(num_inds)
        integer, intent(in) :: ct
        logical :: t_dir

        t_dir = (indices(ct) == indices(ct + step) )
    end function is_direct_pair

    !------------------------------------------------------------------------------------------!    
    
end module LMat_freeze
