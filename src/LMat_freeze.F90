#include "macros.h"

module LMat_freeze
    use IntegralsData, only: nFrozen, UMat, umat_win
    use UMatCache, only: numBasisIndices, UMatInd
    use constants
    use util_mod, only: operator(.div.), custom_findloc
    use OneEInts, only: TMat2D
    use SystemData, only: ECore
    use Parallel_neci
    use ParallelHelper
    implicit none

    private
    public :: freeze_lmat, t_freeze, map_indices, init_freeze_buffers, &
              flush_freeze_buffers, add_core_en

    !> Parameter for number of double excitations a single LMat entry can contribute to
    !! (counting permutations)
    integer, parameter :: num_ex = 4
    integer, parameter :: num_inds = 6
    ! Step corresponding to conjugation of an index
    integer, parameter :: step = (num_inds / 2)

    ! Local storage for ECore and TMat
    real(dp) :: ECore_local, ECore_tot
    HElement_t(dp), allocatable :: TMat_local(:, :), TMat_total(:, :)

    ! Direct orbital excitation = pair of indices in an index array which is apart by step

contains

    !> Initialize the local storage for the diagonal and one-electron terms (two-electron
    !! terms are shared memory). This is called before reading in the 6-index integrals
    subroutine init_freeze_buffers()
        ! Setup the buffers for accumulating the correction per processor

        if(nFrozen > 0) then
            ECore_local = 0.0_dp
            if(allocated(TMat_local)) deallocate(TMat_local)
            allocate(TMat_local(size(TMat2D, dim=1), size(TMat2D, dim=2)))
            if(allocated(TMat_total)) deallocate(TMat_total)
            allocate(TMat_total(size(TMat2D, dim=1), size(TMat2D, dim=2)))
            TMat_local = 0.0_dp
            TMat_total = 0.0_dp

        end if
    end subroutine init_freeze_buffers

    !------------------------------------------------------------------------------------------!

    !> Sum the locally accumulated corrections to the diagonal and one-body terms up and
    !! add them to the global terms, then deallocate temporaries. This is called
    !! after reading in the 6-index integrals
    subroutine flush_freeze_buffers()
        integer(MPIArg) :: tmat_size
        integer(MPIArg) :: ierr

        if(nFrozen > 0) then

            if(.not. allocated(TMat_local) .or. .not. allocated(TMat_total)) &
                call stop_all("flush_freeze_buffers", &
                              "Buffers for freezing three-body interaction not allocated")

            ! Sum the core energy from each proc on this node
            call MPI_Allreduce(ECore_local, ECore_tot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                               mpi_comm_intra, ierr)
            ECore = ECore + ECore_tot

            tmat_size = size(TMat_local)
            ! Sum the 1-body terms from each proc on this node
            call MPI_Allreduce(TMat_local, TMat_total, tmat_size, MPI_DOUBLE_PRECISION, MPI_SUM, &
                               mpi_comm_intra, ierr)
            TMat2D = TMat2D+TMat_total

            ! Deallocate the temporary for the 1-body integrals
            deallocate(TMat_local)
            deallocate(TMat_total)

            ! The 2-body terms are shared memory, no operation is required at this point
        end if
    end subroutine flush_freeze_buffers

    !------------------------------------------------------------------------------------------!

    !> Checks if an entry is zeroed due to frozen orbitals being included
    !> @param[in] indices  array of size 6 containing the indices of the entry in question in the frozen orbital numbering
    !> @return t_freeze  true if the entry is zeroed
    pure function t_freeze(indices)
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
        real(dp) :: prefactor, delta
        integer(MPIArg) :: ierr

        if(t_freeze(indices)) then
            ! Count the number of different indices appearing
            counts = count_frozen_inds(indices)
            ! How many pairs of frozen orbitals do we have
            select case(counts)
            ! 1 => double excitation
            case(1)
                idx = frozen_double_entry(indices, prefactor)
                ! There are up to eight potential double excitations with four given indices
                ! to which the LMat entry contributes - this is because LMat is hermitian
                ! but UMat is not
                do ct = 1, num_ex
                    ! If the index is assigned, there is a duplicate frozen orb
                    if(idx(ct) > 0) then
                        ! Absorb the matrix element into UMat (with the corresponding prefactor
                        ! according to permutation and number of possible spin terms
                        ! UMat(idx(ct)) = UMat(idx(ct)) + prefactor*matel
                        delta = prefactor * matel
                        ! Use a remote update of UMat
                        call MPI_Accumulate(delta, 1, MPI_DOUBLE_PRECISION, 0, &
                                            idx(ct) - 1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, umat_win, ierr)
                    endif
                end do
                ! 2 => single excitation
            case(2)
                idx(1:2) = frozen_single_entry(indices, prefactor)

                if(idx(1) > 0) then
                    ! Absorb the matrix element into TMat
                    ! TMat2D is indexed with spin orbs
                    call add_to_tmat(spatToSpinAlpha(idx(1)), spatToSpinAlpha(idx(2)))
                    call add_to_tmat(spatToSpinBeta(idx(1)), spatToSpinBeta(idx(2)))
                endif
                ! 3 => diagonal element
            case(3)
                call frozen_diagonal_entry(indices, prefactor)
                ECore_local = ECore_local + prefactor * matel
            end select
            ! Zero the matrix element for further usage (i.e. will not turn up anymore)
            matel = 0.0_dp
        endif

    contains

        subroutine add_to_tmat(ind1, ind2)
            integer(int64), intent(in) :: ind1, ind2

            TMat_local(ind1, ind2) = TMat_local(ind1, ind2) + prefactor * matel
            ! TMat is hermitian, as is the 3-body term
            if(ind1 /= ind2) &
                TMat_local(ind2, ind1) = TMat_local(ind2, ind1) + prefactor * matel
        end subroutine add_to_tmat

    end subroutine add_core_en

    !------------------------------------------------------------------------------------------!

    !> Count how many frozen indices there are, determine if the entry has to be added
    !! to a lower-order term and if yes, which one.
    !> @param[in] indices  array of size 6 with the orbital indices of the entry
    !> @return level  indicates to which lower order term this entry has to be added
    !!                and therefore, which subroutine has to be called in the following
    !!                0 - none
    !!                1 - one frozen index pair -> two-body term
    !!                2 - two frozen index pairs -> one-body term
    !!                3 - three frozen index pairs -> diagonal term
    pure function count_frozen_inds(indices) result(level)
        integer(int64), intent(in) :: indices(num_inds)
        integer :: level
        integer :: counts

        integer :: ct
        level = 0
        ! If there is an unpaired frozen orbital, the entry is discarded, it does not
        ! contribute to any contraction
        do ct = 1, num_inds
            if(count(indices(ct) == indices) == 1 .and. indices(ct) < 1) return
        end do
        ! Else, it contributes to a contraction according to the number of frozen indices
        counts = count(indices < 1)
        ! An even number of frozen orbital indices is required
        if(modulo(counts, 2) /= 0) return

        level = counts.div.2

    end function count_frozen_inds

    !------------------------------------------------------------------------------------------!

    !> Get the index of the UMat entry to which the LMat entry with given indices shall be
    !! added if it is frozen.
    !> @param[in] indices  array of lenght 6 with the orbital indices of the LMat entry
    !> @param[ou] t_par  flag indicating if the matrix element enters UMat with a -1
    !> @return index  index of the UMat entry to add the LMat entry to, 0 if entry is not frozen
    function frozen_double_entry(indices, prefactor) result(idx)
        integer(int64), intent(in) :: indices(num_inds)
        real(dp), intent(out) :: prefactor
        integer(int64) :: idx(num_ex)
        integer :: f_one, f_two, ct, counts, i, rs, marks(3)
        logical :: unfrozen(num_inds)
        integer(int64) :: uf_idx(4)

        idx = 0
        unfrozen = indices > 0
        ! Mark where the unfrozen indices are
        marks = 0

        ! Identify the direct unfrozen pair
        ! Look for the direct pair of unfrozen indices (these shall stay direct)
        do ct = 1, step
            if(unfrozen(ct) .and. unfrozen(ct + step)) then
                uf_idx(1) = indices(ct)
                uf_idx(3) = indices(ct + step)
                ! Mark the extracted indices so they are not taken again
                marks(1) = ct
                marks(2) = ct + step
                exit
            end if
        end do
        ! Now, find the other two (direct or exchange) unfrozen indices
        rs = 2
        do ct = 1, num_inds
            if(unfrozen(ct) .and. .not. any(ct == marks)) then
                uf_idx(rs) = indices(ct)
                marks(3) = ct
                rs = rs + 2
            endif
        end do

        ! Now, get the required UMat indices, including all transposed index pairs
        idx = permute_umat_inds(int(uf_idx(1)), int(uf_idx(2)), int(uf_idx(3)), int(uf_idx(4)))
        ! Check the permutation and possible factor of two due to spin
        ! First, get the two frozen orbitals
        f_one = custom_findloc(unfrozen, .false., back=.false.)
        f_two = custom_findloc(unfrozen, .false., back=.true.)

        ! There are two spin configurations in a close-shell frozen scenario for terms
        ! with a direct frozen orbital (i.e. no permutation required)
        ! All others enter with a prefactor of -1
        if(is_direct(f_one, f_two)) then
            prefactor = 2.0_dp
        else
            prefactor = -1.0_dp
            ! If this is a diagonal term (only three different indices), permutational symmetry
            ! gives an extra factor of 2 for the exchange term. If only two different
            ! indices appear, also the direct term gets this prefactor
            counts = 1
            ! Count the number of different indices
            do ct = 2, num_inds
                if(.not. any(indices(ct) == indices(1:ct - 1))) counts = counts + 1
            end do
            if(counts == 2) then
                prefactor = 2.0_dp * prefactor
            elseif(counts == 3) then
                ! If there are three different indices, only add the extra factor if there is no direct unfrozen pair
                do ct = 1, step
                    if(is_repeated_pair(indices, ct) .and. unfrozen(ct)) return
                end do
                prefactor = 2.0_dp * prefactor
            end if
        end if

    end function frozen_double_entry

    !------------------------------------------------------------------------------------------!

    !> Returns the UMatInd values of all possible permutations of the input indices
    !> @param[in] a,b,c,d  orbital indices of a two-body element
    !> @return inds  array of size 4 containing the indices of UMat corresponding to these
    !!               four orbitals (in this order) and their hermitian conjugates
    !!               entries of 0 indicate that no position in UMat has to be addressed
    function permute_umat_inds(a, b, c, d) result(inds)
        integer, intent(in) :: a, b, c, d
        integer(int64) :: inds(num_ex)

        integer :: ct

        inds = 0
        ! These are adjoint to each other, they are the same in LMat, but not UMat
        inds(1) = UMatInd(a, b, c, d)
        inds(2) = UMatInd(c, b, a, d)
        inds(3) = UMatInd(a, d, c, b)
        inds(4) = UMatInd(c, d, a, b)

        ! All terms are only counted if they are different from a previously occuring index
        ! This approach might be less performant than competing ways of doing this check,
        ! but it is very clear and less prone to bugs. If this proves to be a performance
        ! bottleneck, it can be changed to a faster implementation later on
        do ct = 2, num_ex
            ! If the index already appeared, delete it
            if(any(inds(ct) == inds(1:ct - 1))) inds(ct) = 0
        end do
    end function permute_umat_inds

    !------------------------------------------------------------------------------------------!

    !> For a given set of 6 orbital indices with two pairs of frozen orbitals,
    !! returns the orbital indices of the corresponding single excitation and the prefactor
    !! for the matrix element due to spin
    !> @param[in] indices  array of size 6 with orbital indices, four of which have to be
    !!                     repeated frozen indices
    !> @param[out] prefactor  on return, the prefactor of the matrix element when added to
    !!                        the one-body terms
    !> @return orbs  the two non-frozen orbitals
    function frozen_single_entry(indices, prefactor) result(orbs)
        integer(int64), intent(in) :: indices(num_inds)
        real(dp), intent(out) :: prefactor

        integer :: orbs(2)
        integer :: f_orbs(num_inds - 2), ct
        integer :: uf_one, uf_two, directs
        logical :: unfrozen(num_inds)

        orbs = 0
        unfrozen = indices > 0
        f_orbs = int(pack(indices,.not. unfrozen))
        ! For the single entry, there are only two unfrozen orbs, find these and return them
        orbs = int(pack(indices, unfrozen))

        ! Again, two things to check: The number of direct orbitals and the permutation

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
            directs = 0
            ! Count the number of direct excitations
            do ct = 1, step
                ! Each direct repeated orbital doubles the number of spin configs
                ! (the case of four identical orbs is excluded)
                if(is_repeated_pair(indices, ct) .or. (unfrozen(ct) .and. unfrozen(ct + step))) then
                    directs = directs + 1
                endif
            end do
            select case(directs)
                ! If there are none, there is no freedom for spin choice, but the LMat entry
                ! appears twice in the matrix element evaluation due to symmetry if the
                ! unfrozen entry is repeated
            case(0)
                if(orbs(1) == orbs(2)) prefactor = 2.0_dp
                ! If there is one, two possible spin configs contribute, and we have odd parity
            case(1)
                prefactor = -2.0_dp
            case(3)
                ! If there are three (just two is not possible), there are two spin degrees of
                ! freedom -> factor of 4
                prefactor = 4.0_dp
            end select
        else
            ! b) All frozen orbs are the same -> only one spin config is allowed, only parity of
            ! the permutation counts
            ! Here, only a direct exctiation has even parity (same rule as above)

            ! Get the posisitons of the unfrozen orbs
            uf_one = custom_findloc(unfrozen, .true., back=.false.)
            uf_two = custom_findloc(unfrozen, .true., back=.true.)

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
        t_quad = .false.
        do ct = 1, 3
            if(count(indices(ct) == indices) > 4) then
                prefactor = 0.0_dp
                return
            end if
            if(count(indices(ct) == indices) > 3) t_quad = .true.
        end do

        ! Both relate to the number of different direct excits, so count these
        directs = 0
        do ct = 1, step
            if(is_repeated_pair(indices, ct)) directs = directs + 1
        end do

        ! If there is a quadruple index, the spin of it is fixed (both have to be occupied), then,
        ! the prefactor is always +/- 2, depending on if there is more than one direct pair
        if(t_quad) then
            ! The only options with a quadruple index are
            ! a) one direct excitation => -2
            ! b) three direct exctiations => +2            
            prefactor = merge(-2.0_dp, 2.0_dp, directs == 1)
        else
            ! In the other case, there are three relevant cases:
            ! a) All excitations are direct => prefactor 8 = 2**3 (all spins are free)
            ! b) One excitation is direct => prefactor -4 (from the spin freedom of the direct orb + overall spin)
            ! c) No excitation is direct => prefactor 4
            !    the factor of 4 comes from two sources: a factor of 2 from the overall spin
            !    and another factor of two from permutational symmetry: both indirect terms
            !    are identical in this case, so this entry has to be counted as both
            if(directs == 0) then
                prefactor = 4.0_dp
            elseif(directs == 1) then
                prefactor = -4.0_dp
            else ! directs == 3
                prefactor = 8.0_dp
            endif
        endif

    end subroutine frozen_diagonal_entry

    !------------------------------------------------------------------------------------------!

    !> For two positions of indices in a 6-index set, return if these are a direct pair
    !> @param[in] one, two  integers between 1 and 6
    !> @return t_dir  true if the two integers correspond to positions that are direct, i.e.
    !!                for which the LMat entry is symmetric under exchange of the values
    pure function is_direct(one, two) result(t_dir)
        integer, intent(in) :: one, two
        logical :: t_dir

        t_dir = (one + step == two)
    end function is_direct

    !------------------------------------------------------------------------------------------!

    !> For a set of 6 orbital indices, return whether the pair at a given position is repeated
    !> @param[in] indices  array of size 6 containing a set of indices indexing an LMat entry
    !> @param[in] ct  integer between 1 and 3 labeling a position in the index set
    !> @return t_dir  true if the two indices at the position ct (i.e. ct and ct+step in indices)
    !!                are the same
    pure function is_repeated_pair(indices, ct) result(t_dir)
        integer(int64), intent(in) :: indices(num_inds)
        integer, intent(in) :: ct
        logical :: t_dir

        t_dir = (indices(ct) == indices(ct + step))
    end function is_repeated_pair

    !------------------------------------------------------------------------------------------!

end module LMat_freeze
