#include "macros.h"
module pchb_excitgen
    use constants
    use SystemData, only: nel, nBasis, G1, tStoreSpinOrbs, AB_elec_pairs, par_elec_pairs, &
                          t_pchb_weighted_singles
    use bit_rep_data, only: NIfTot
    use dSFMT_interface, only: genrand_real2_dSFMT
    use get_excit, only: make_double, exciteIlut
    use excit_gens_int_weighted, only: pick_biased_elecs, weighted_single_excit_wrapper, &
                                       pgen_single_4ind
    use FciMCData, only: pSingles, excit_gen_store_type, nInvalidExcits, nValidExcits, &
                         pParallel, pDoubles
    use excitation_types, only: DoubleExc_t
    use sltcnd_mod, only: sltcnd_excit
    use procedure_pointers, only: generate_single_excit
    use UMatCache, only: gtID, numBasisIndices
    use aliasSampling, only: aliasSamplerArray_t
    use util_mod, only: fuseIndex, linearIndex, intswap, getSpinIndex, near_zero
    use GenRandSymExcitNUMod, only: uniform_single_excit_wrapper, calc_pgen_symrandexcit2
    use SymExcitDataMod, only: pDoubNew, scratchSize
    use sym_general_mod, only: IsSymAllowedExcitMat
    implicit none

    ! there are three pchb_samplers:
    ! 1 - same-spin case
    ! 2 - opp spin case without exchange
    ! 3 - opp spin case with exchange
    integer, parameter :: SAME_SPIN = 1, OPP_SPIN_NO_EXCH = 2, OPP_SPIN_EXCH = 3
    type(aliasSamplerArray_t) :: pchb_samplers(3)
    real(dp), allocatable :: pExch(:)
    integer, allocatable :: tgtOrbs(:, :)

contains

    ! this is the interface routine: for singles, use the uniform excitgen
    ! for doubles, the precomputed heat-bath weights
    subroutine gen_rand_excit_pchb(nI, ilutI, nJ, ilutJ, exFlag, ic, ex, tpar, &
                                   pgen, helgen, store, part_type)
        implicit none
        ! The interface is common to all excitation generators, see proc_ptrs.F90
        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ic, ex(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:NIfTot)
        logical, intent(out) :: tpar
        real(dp), intent(out) :: pGen
        HElement_t(dp), intent(out) :: HElGen
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: part_type

        unused_var(exFlag); unused_var(part_type)

        helgen = 0.0_dp

        if (genrand_real2_dSFMT() < pSingles) then
            ic = 1
            ! defaults to uniform singles, but can be set to other excitgens
            call generate_single_excit(nI, ilutI, nJ, ilutJ, ex, tpar, store, pgen)
        else
            ic = 2
            ! use precomputed weights to generate doubles
            call generate_double_pchb(nI, ilutI, nJ, ilutJ, ex, tpar, pgen)
            pgen = pgen * (1.0 - pSingles)

            if (IsNullDet(nJ)) then
                nInvalidExcits = nInvalidExcits + 1
            else
                nValidExcits = nValidExcits + 1
            end if
        end if

    end subroutine gen_rand_excit_pchb

    !------------------------------------------------------------------------------------------!

    !> given the initial determinant (both as nI and ilut), create a random double
    !! excitation using the hamiltonian matrix elements as weights
    !> @param[in] nI  determinant to excite from
    !> @param[in] elec_map  map to translate electron picks to orbitals
    !> @param[in] ilut  determinant to excite from in ilut format
    !> @param[out] nJ  on return, excited determinant
    !> @param[out] excitMat  on return, excitation matrix nI -> nJ
    !> @param[out] tParity  on return, the parity of the excitation nI -> nJ
    !> @param[out] pGen  on return, the probability of generating the excitation nI -> nJ
    subroutine generate_double_pchb(nI, ilutI, nJ, ilutJ, ex, tpar, pgen)
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel)
        integer(n_int), intent(out) :: ilutJ(0:NIfTot)
        integer, intent(out) :: ex(2, maxExcit)
        real(dp), intent(out) :: pGen
        logical, intent(out) :: tpar

        integer :: elecs(2), src(2), sym_prod, ispn, sum_ml, ij
        integer :: orbs(2), ab
        real(dp) :: pGenHoles
        logical :: invalid
        integer :: spin(2), samplerIndex

        ! first, pick two random elecs
        call pick_biased_elecs(nI, elecs, src, sym_prod, ispn, sum_ml, pGen)
        if (src(1) > src(2)) call intswap(src(1), src(2))

        invalid = .false.
        ! use the sampler for this electron pair -> order of src electrons does not matter
        ij = fuseIndex(gtID(src(1)), gtID(src(2)))
        ! the spin of the electrons: 0 - alpha, 1 - beta
        spin = getSpinIndex(src)
        ! determine type of spin-excitation: same-spin, opp spin w exchange, opp spin w/o exchange
        if (spin(1) == spin(2)) then
            ! same spin
            samplerIndex = SAME_SPIN
        else
            ! else, pick exchange with...some ij-spin bias
            if (genrand_real2_dSFMT() < pExch(ij)) then
                samplerIndex = OPP_SPIN_EXCH
                ! adjust pgen
                pGen = pGen * pExch(ij)
                ! the spins of the target are the opposite of the source spins
                call intswap(spin(1), spin(2))
            else
                samplerIndex = OPP_SPIN_NO_EXCH
                ! adjust pgen
                pGen = pGen * (1.0_dp - pExch(ij))
            end if
        end if
        ! get a pair of orbitals using the precomputed weights
        call pchb_samplers(samplerIndex)%aSample(ij, ab, pGenHoles)
        ! split the index ab (using a table containing mapping ab -> (a,b))
        orbs = tgtOrbs(:, ab)
        ! convert orbs to spin-orbs with the same spin
        orbs = 2 * orbs - spin

        ! check if the picked orbs are a valid choice - if they are the same, match one
        ! occupied orbital or are zero (maybe because there are no allowed picks for
        ! the given source) abort
        invalid = (any(orbs == 0) .or. any(orbs(1) == nI) &
                   .or. any(orbs(2) == nI)) .or. (orbs(1) == orbs(2))
        ! unfortunately, there is a super-rare case when, due to floating point error,
        ! an excitation with pGen=0 is created. Those are invalid, too
        if (near_zero(pGenHoles)) then
            invalid = .true.
            ! Yes, print. Those events are signficant enough to be always noted in the output
            print *, "WARNING: Generated excitation with probability of 0"
        end if

        pGen = pGen * pGenHoles
        if (invalid) then
            ! if 0 is returned, there are no excitations for the chosen elecs
            ! -> return nulldet
            nJ = 0
            ilutJ = 0_n_int
            ex(2, 1:2) = orbs
            ex(1, 1:2) = src
        else
            ! else, construct the det from the chosen orbs/elecs

            call make_double(nI, nJ, elecs(1), elecs(2), orbs(1), orbs(2), ex, tpar)

            ilutJ = exciteIlut(ilutI, src, orbs)
        end if

    end subroutine generate_double_pchb

    !------------------------------------------------------------------------------------------!

    !> Calculate the probability of generating a given excitation with the pchb excitgen
    !> @param[in] nI  determinant to start from
    !> @param[in] ex  2x2 excitation matrix
    !> @param[in] ic  excitation level
    !> @param[in] ClassCount2  symmetry information of the determinant
    !> @param[in] ClassCountUnocc2  symmetry information of the virtual orbitals
    !> @return pGen  probability of drawing this excitation with the pchb excitgen
    function calc_pgen_pchb(nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2) result(pgen)
        implicit none
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(in) :: ex(2, 2), ic
        integer, intent(in) :: ClassCount2(ScratchSize), ClassCountUnocc2(ScratchSize)
        real(dp) :: pgen

        if (ic == 1) then
            ! single excitations are the job of the uniform excitgen
            if (t_pchb_weighted_singles) then
                pgen = pSingles * pgen_single_4ind(nI, ilutI, ex(1, 1), ex(2, 1))
            else
                call calc_pgen_symrandexcit2(nI, ex, ic, ClassCount2, ClassCountUnocc2, pDoubles, pGen)
            end if
        else if (ic == 2) then
            pGen = (1.0 - pSingles) * calc_double_pgen_pchb(ex)
        else
            pgen = 0.0_dp
        end if

    end function calc_pgen_pchb

    !------------------------------------------------------------------------------------------!

    !> Calculate the probability of drawing a given double excitation ex
    !> @param[in] ex  2x2 excitation matrix
    !> @return pgen  probability of generating this double with the pchb double excitgen
    function calc_double_pgen_pchb(ex) result(pgen)
        implicit none
        integer, intent(in) :: ex(2, 2)
        real(dp) :: pgen
        integer :: ab, ij, nex(2, 2), samplerIndex

        ! spatial orbitals of the excitation
        nex = gtID(ex)
        ij = fuseIndex(nex(1, 1), nex(1, 2))
        ! the probability of picking the two electrons: they are chosen uniformly
        ! check which sampler was used
        if (is_beta(ex(1, 1)) .eqv. is_beta(ex(1, 2))) then
            pgen = pParallel / par_elec_pairs
            ! same-spin case
            samplerIndex = SAME_SPIN
        else
            pgen = (1.0_dp - pParallel) / AB_elec_pairs
            ! excitations without spin-exchange OR to the same spatial orb
            if ((is_beta(ex(1, 1)) .eqv. is_beta(ex(2, 1))) .or. (nex(2, 1) == nex(2, 2))) then
                ! opp spin case without exchange
                samplerIndex = OPP_SPIN_NO_EXCH
                pgen = pgen * (1 - pExch(ij))
            else
                ! opp spin case with exchange
                samplerIndex = OPP_SPIN_EXCH
                pgen = pgen * pExch(ij)
            end if
        end if

        ! look up the probability for this excitation in the sampler
        ab = fuseIndex(nex(2, 1), nex(2, 2))
        pgen = pgen * pchb_samplers(samplerIndex)%aGetProb(ij, ab)

    end function calc_double_pgen_pchb

    !------------------------------------------------------------------------------------------!

    !> initialize the pchb excitation generator
    !! this does two things:
    !! 1. setup the lookup table for the mapping ab -> (a,b)
    !! 2. setup the alias table for picking ab given ij with probability ~<ij|H|ab>
    subroutine init_pchb_excitgen(projEDet)
        implicit none
        integer, intent(in) :: projEDet(:)
        integer :: ab, a, b, abMax
        integer :: aerr, nBI
        integer(int64) :: memCost
        integer :: samplerIndex

        write(iout, *) "Allocating PCHB excitation generator objects"
        ! total memory cost
        memCost = 0_int64
        ! number of spatial orbs
        nBI = numBasisIndices(nBasis)
        ! initialize the mapping ab -> (a,b)
        abMax = fuseIndex(nBI, nBI)
        allocate(tgtOrbs(2, 0:abMax), stat=aerr)
        do a = 1, nBI
            do b = 1, a
                ab = fuseIndex(a, b)
                tgtOrbs(1, ab) = b
                tgtOrbs(2, ab) = a
            end do
        end do

        ! enable catching exceptions
        tgtOrbs(:, 0) = 0

        ! setup the alias table
        call setup_pchb_sampler(projEDet)

        write(iout, *) "Finished excitation generator initialization"
        ! this is some bias used internally by CreateSingleExcit - not used here
        pDoubNew = 0.0

        ! Set the single excitation generator
        if (t_pchb_weighted_singles) then
            generate_single_excit => weighted_single_excit_wrapper
        else
            generate_single_excit => uniform_single_excit_wrapper
        end if

    contains

        subroutine setup_pchb_sampler(projEDet)
            implicit none
            integer, intent(in) :: projEDet(:)
            integer :: i, j, iSampler
            integer :: ij, ijMax
            integer :: ex(2, 2)
            real(dp), allocatable :: w(:)
            real(dp), allocatable :: pNoExch(:)
            ! number of possible source orbital pairs
            ijMax = fuseIndex(nBI, nBI)
            ! allocate the bias for picking an exchange excitation
            allocate(pExch(ijMax), stat=aerr, source=0.0_dp)
            ! temporary storage for the unnormalized prob of not picking an exchange excitation
            allocate(pNoExch(ijMax), stat=aerr, source=1.0_dp)

            memCost = memCost + abMax * ijMax * 24 * 3
            write(iout, *) "Excitation generator requires", real(memCost, dp) / 2.0_dp**30, "GB of memory"
            write(iout, *) "Generating samplers for PCHB excitation generator"
            ! weights per pair
            allocate(w(abMax), stat=aerr)
            ! initialize the three samplers
            do samplerIndex = 1, 3
                ! allocate: all samplers have the same size
                call pchb_samplers(samplerIndex)%setupSamplerArray(int(ijMax, int64), int(abMax, int64))
                do i = 1, nBI
                    ! map i to alpha spin (arbitrary choice)
                    ex(1, 1) = 2 * i
                    ! as we order a,b, we can assume j <= i
                    do j = 1, i
                        w(:) = 0.0_dp
                        ! for samplerIndex == 1, j is alpha, else, j is beta
                        ex(1, 2) = map_orb(j, (/1/))
                        ! for each (i,j), get all matrix elements <ij|H|ab> and use them as
                        ! weights to prepare the sampler
                        do a = 1, nBI
                            ! a is alpha for same-spin (1) and opp spin w/o exchange (2)
                            ex(2, 2) = map_orb(a, (/1, 2/))
                            do b = 1, a
                                ! exception: for sampler 3, a!=b
                                if (samplerIndex == OPP_SPIN_EXCH .and. a == b) cycle
                                ab = fuseIndex(a, b)
                                ! ex(2,:) is in ascending order
                                ! b is alpha for sampe-spin (1) and opp spin w exchange (3)
                                ex(2, 1) = map_orb(b, (/1, 3/))
                                ! use the actual matrix elements as weights
                                w(ab) = abs(sltcnd_excit(projEDet, DoubleExc_t(ex), .false.))
                            end do
                        end do
                        ij = fuseIndex(i, j)
                        call pchb_samplers(samplerIndex)%setupEntry(ij, w)
                        if (samplerIndex == OPP_SPIN_EXCH) pExch(ij) = sum(w)
                        if (samplerIndex == OPP_SPIN_NO_EXCH) pNoExch(ij) = sum(w)
                    end do
                end do
            end do

            ! normalize the exchange bias (where normalizable)
            where (.not. near_zero(pExch + pNoExch)) pExch = pExch / (pExch + pNoExch)

            deallocate(w)
            deallocate(pNoExch)
        end subroutine setup_pchb_sampler

        function map_orb(orb, alphaSamplers) result(sorb)
            ! map spatial orbital to the spin orbital matching the current samplerIndex
            ! Input: orb - spatial orbital to be mapped
            !        alphaSamplers - list of samplerIndex values for which the mapping shall be to alpha
            ! Output: sorb - corresponding spin orbital
            integer, intent(in) :: orb
            integer, intent(in) :: alphaSamplers(:)
            integer :: sorb

            if (any(samplerIndex == alphaSamplers)) then
                sorb = 2 * orb
            else
                sorb = 2 * orb - 1
            end if
        end function map_orb

    end subroutine init_pchb_excitgen

    !------------------------------------------------------------------------------------------!

    !> deallocate the sampler and the mapping ab -> (a,b)
    subroutine finalize_pchb_excitgen()
        implicit none
        integer :: samplerIndex

        do samplerIndex = 1, 3
            call pchb_samplers(samplerIndex)%samplerArrayDestructor()
        end do
        deallocate(tgtOrbs)
        deallocate(pExch)

    end subroutine finalize_pchb_excitgen

end module pchb_excitgen
