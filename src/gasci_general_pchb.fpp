#include "macros.h"
#:include "macros.fpph"
#:include "algorithms.fpph"

module gasci_general_pchb
    use constants, only: n_int, dp, int64, maxExcit, iout
    use gasci_general, only: gen_exc_single
    use util_mod, only: fuseIndex, getSpinIndex, near_zero, intswap, operator(.div.)
    use dSFMT_interface, only: genrand_real2_dSFMT
    use get_excit, only: make_double, exciteIlut
    use SymExcitDataMod, only: pDoubNew, ScratchSize
    use excitation_types, only: SingleExc_t, DoubleExc_t
    use sltcnd_mod, only: sltcnd_excit
    use procedure_pointers, only: generate_single_excit_t
    use aliasSampling, only: AliasSampler_3D_t
    use UMatCache, only: gtID, numBasisIndices
    use FciMCData, only: pSingles, excit_gen_store_type, nInvalidExcits, nValidExcits, &
                         pParallel, projEDet
    use excit_gens_int_weighted, only: pick_biased_elecs

    use SystemData, only: nEl, AB_elec_pairs, par_elec_pairs
    use bit_rep_data, only: NIfTot

    use gasci, only: GASSpec_t
    use gasci_supergroup_index, only: SuperGroupIndexer_t

    implicit none

    private
    public :: gen_GASCI_general_pchb, general_GAS_PCHB, GAS_PCHB_excit_gen_t

    abstract interface
        real(dp) function calc_pgen_t(nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2)
            import :: n_int, dp, NIfTot, nEl, ScratchSize
            integer, intent(in) :: nI(nel)
            integer(n_int), intent(in) :: ilutI(0:NIfTot)
            integer, intent(in) :: ex(2, 2), ic
            integer, intent(in) :: ClassCount2(ScratchSize), ClassCountUnocc2(ScratchSize)
        end function
    end interface

    ! there are three pchb_samplers for each supergroup:
    ! 1 - same-spin case
    ! 2 - opp spin case without exchange
    ! 3 - opp spin case with exchange
    integer, parameter :: SAME_SPIN = 1, OPP_SPIN_NO_EXCH = 2, OPP_SPIN_EXCH = 3
!
    type :: GAS_PCHB_excit_gen_t
        private
        !> The shape is (fused_number_of_double_excitations, 3, n_supergroup)
        type(AliasSampler_3D_t) :: pchb_samplers
        type(SuperGroupIndexer_t) :: indexer
        type(GASSpec_t) :: GASSpec
        real(dp), allocatable :: pExch(:, :)
        integer, allocatable :: tgtOrbs(:, :)

        ! I would prefer to initialize directly as
        !       generate_single => generate_single_GAS
        ! and use
        !       if (present(generate_single)) this%generate_single => generate_single
        ! in the init routines.
        ! Unfortunately this causes a **compile time** segfault in ifort <= 18.x.x
        procedure(generate_single_excit_t), pointer, nopass :: generate_single => null()
        procedure(calc_pgen_t), pointer, nopass :: calc_pgen_single => null()

    contains
        private
        procedure, public :: init => init_pchb_excitgen
        procedure, public :: finalize
        procedure, public :: gen_excit
        procedure, public :: calc_pgen => calc_pgen_pchb

        procedure :: calc_double_pgen_pchb
        procedure :: generate_double
        procedure :: is_allowed
    end type

    type(GAS_PCHB_excit_gen_t) :: general_GAS_PCHB

contains

    !>  @brief
    !>  The disconnected_GAS_PCHB excitation generator subroutine.
    !>
    !>  @details
    !>  This is a wrapper around `disconnected_GAS_PCHB%gen_excit`
    !>  to match the function pointer interface.
    !>  The interface is common to all excitation generators, see proc_ptrs.F90
    subroutine gen_GASCI_general_pchb(nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                        ex_mat, tParity, pGen, hel, store, part_type)
        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ic, ex_mat(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:NifTot)
        real(dp), intent(out) :: pGen
        logical, intent(out) :: tParity
        HElement_t(dp), intent(out) :: hel
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: part_type
        character(*), parameter :: this_routine = 'gen_GASCI_pchb'

        @:unused_var(exFlag, part_type, store)
        @:ASSERT(general_GAS_PCHB%GASSpec%contains_det(nI))

        hel = h_cast(0.0_dp)
        call general_GAS_PCHB%gen_excit(nI, ilutI, nJ, ilutJ, ic, ex_mat, tParity, store, pgen)
    end subroutine gen_GASCI_general_pchb

    !> @brief
    !> Check if a double excitation is allowed.
    logical pure function is_allowed(self, exc, supergroup)
        class(GAS_PCHB_excit_gen_t), intent(in) :: self
        type(DoubleExc_t), intent(in) :: exc
        integer, intent(in) :: supergroup(:)

        integer :: excited_supergroup(size(supergroup))
        integer :: src_spaces(2), tgt_spaces(2), i

        src_spaces = self%GASSpec%get_iGAS(exc%val(1, :))
        tgt_spaces = self%GASSpec%get_iGAS(exc%val(2, :))

        if (all(src_spaces == tgt_spaces) .and. src_spaces(1) == src_spaces(2)) then
            ! All electrons come from the same space and there are no restrictions
            ! regarding recoupling.
            is_allowed = .true.
        else
            ! Ensure that GAS specifications contain supergroup **after** excitation.
            excited_supergroup = supergroup
            excited_supergroup(src_spaces) = excited_supergroup(src_spaces) - 1
            excited_supergroup(tgt_spaces) = excited_supergroup(tgt_spaces) + 1

            is_allowed = self%GASSpec%contains_supergroup(excited_supergroup)
        end if
    end function

!
    !>  @brief
    !>  Initialize the pchb excitation generator
    !>
    !>  @details
    !>  This does two things:
    !>  1. setup the lookup table for the mapping ab -> (a,b)
    !>  2. setup the alias table for picking ab given ij with probability ~<ij|H|ab>
    subroutine init_pchb_excitgen(this, GASSpec, generate_single, calc_pgen_single)
        class(GAS_PCHB_excit_gen_t), intent(out) :: this
        type(GASSpec_t), intent(in) :: GASSpec
        procedure(generate_single_excit_t), optional :: generate_single
        procedure(calc_pgen_t), optional :: calc_pgen_single

        integer :: ab, a, b, abMax
        integer :: aerr, nBI
        integer(int64) :: memCost
        integer :: samplerIndex

        this%indexer = SuperGroupIndexer_t(GASSpec)
        this%GASSpec = GASSpec

        write(iout, *) "Allocating PCHB excitation generator objects"
        ! total memory cost
        memCost = 0_int64
        ! number of spatial orbs
        nBI = numBasisIndices(this%GASSpec%n_spin_orbs())
        ! initialize the mapping ab -> (a,b)
        abMax = fuseIndex(nBI, nBI)
        allocate(this%tgtOrbs(2, 0 : abMax))
        do a = 1, nBI
            do b = 1, a
                ab = fuseIndex(a, b)
                this%tgtOrbs(1, ab) = b
                this%tgtOrbs(2, ab) = a
            end do
        end do

        ! enable catching exceptions
        this%tgtOrbs(:, 0) = 0

        ! setup the alias table
        call setup_pchb_sampler()

        write(iout, *) "Finished excitation generator initialization"
        ! this is some bias used internally by CreateSingleExcit - not used here
        pDoubNew = 0.0

        if (present(generate_single)) then
            this%generate_single => generate_single
        else
            this%generate_single => generate_single_GAS
        end if
        if (present(calc_pgen_single)) then
            this%calc_pgen_single => calc_pgen_single
        else
            this%calc_pgen_single => calc_pgen_single_todo
        end if
    contains

        subroutine setup_pchb_sampler()
            integer :: i, j, iSampler
            integer :: ij, ijMax
            integer :: ex(2, 2)
            real(dp), allocatable :: w(:)
            real(dp), allocatable :: pNoExch(:)
            integer, allocatable :: supergroups(:, :)
            integer :: i_sg
            ! possible supergroups
            supergroups = this%indexer%get_supergroups()
            ! number of possible source orbital pairs
            ijMax = fuseIndex(nBI, nBI)
            ! allocate the bias for picking an exchange excitation
            allocate(this%pExch(ijMax, size(supergroups, 2)), source=0.0_dp)
            ! temporary storage for the unnormalized prob of not picking an exchange excitation

            memCost = size(supergroups, 2) * (memCost + abMax * ijMax * 24 * 3)

            call this%pchb_samplers%shared_alloc([ijMax, 3, size(supergroups, 2)], abMax)
            write(iout, *) "Excitation generator requires", real(memCost, dp) / 2.0_dp**30, "GB of memory"
            write(iout, *) "Generating samplers for PCHB excitation generator"
            ! weights per pair
            allocate(w(abMax))
            ! initialize the three samplers
            do i_sg = 1, size(supergroups, 2)
                if (mod(i_sg, 100) == 0) write(iout, *) 'Still generating the samplers'
                pNoExch = 1.0_dp - this%pExch(:, i_sg)
                do samplerIndex = 1, 3
                    ! allocate: all samplers have the same size
                    do i = 1, nBI
                        ! map i to alpha spin (arbitrary choice)
                        ex(1, 1) = 2 * i
                        ! as we order a,b, we can assume j <= i
                        do j = 1, i
                            w(:) = 0.0_dp
                            ! for samplerIndex == 1, j is alpha, else, j is beta
                            ex(1, 2) = map_orb(j, [SAME_SPIN])
                            ! for each (i,j), get all matrix elements <ij|H|ab> and use them as
                            ! weights to prepare the sampler
                            do a = 1, nBI
                                ! a is alpha for same-spin (1) and opp spin w/o exchange (2)
                                ex(2, 2) = map_orb(a, [SAME_SPIN, OPP_SPIN_NO_EXCH])
                                do b = 1, a
                                    ! exception: for sampler 3, a!=b
                                    if (samplerIndex == OPP_SPIN_EXCH .and. a == b) cycle
                                    ab = fuseIndex(a, b)
                                    ! ex(2,:) is in ascending order
                                    ! b is alpha for sampe-spin (1) and opp spin w exchange (3)
                                    ex(2, 1) = map_orb(b, [SAME_SPIN, OPP_SPIN_EXCH])
                                    ! use the actual matrix elements as weights
                                    if (this%is_allowed(DoubleExc_t(ex), supergroups(:, i_sg))) then
                                        w(ab) = abs(sltcnd_excit(projEDet(:, 1), DoubleExc_t(ex), .false.))
                                    else
                                        w(ab) = 0._dp
                                    end if
                                end do
                            end do
                            ij = fuseIndex(i, j)
                            call this%pchb_samplers%setup_entry(ij, samplerIndex, i_sg, w)
                            if (samplerIndex == OPP_SPIN_EXCH) this%pExch(ij, i_sg) = sum(w)
                            if (samplerIndex == OPP_SPIN_NO_EXCH) pNoExch(ij) = sum(w)
                        end do
                    end do
                end do
                ! normalize the exchange bias (where normalizable)
                where (near_zero(this%pExch(:, i_sg) + pNoExch))
                    this%pExch(:, i_sg) = 0._dp
                else where
                    this%pExch(:, i_sg) = this%pExch(:, i_sg) / (this%pExch(:, i_sg) + pNoExch)
                end where
            end do

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

    !>  @brief
    !>  Given the initial determinant (both as nI and ilut), create a random double
    !>  excitation using the hamiltonian matrix elements as weights
    !>
    !> @param[in] nI  determinant to excite from
    !> @param[in] elec_map  map to translate electron picks to orbitals
    !> @param[in] ilut  determinant to excite from in ilut format
    !> @param[out] nJ  on return, excited determinant
    !> @param[out] excitMat  on return, excitation matrix nI -> nJ
    !> @param[out] tParity  on return, the parity of the excitation nI -> nJ
    !> @param[out] pGen  on return, the probability of generating the excitation nI -> nJ
    subroutine generate_double(this, nI, ilutI, nJ, ilutJ, ex, tpar, pgen)
        class(GAS_PCHB_excit_gen_t), intent(in) :: this
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
        integer(int64) :: i_sg
        debug_function_name("generate_double")


        ! first, pick two random elecs
        call pick_biased_elecs(nI, elecs, src, sym_prod, ispn, sum_ml, pGen)
        if (src(1) > src(2)) call intswap(src(1), src(2))

        i_sg = this%indexer%idx_nI(nI)

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
            if (genrand_real2_dSFMT() < this%pExch(ij, i_sg)) then
                samplerIndex = OPP_SPIN_EXCH
                ! adjust pgen
                pGen = pGen * this%pExch(ij, i_sg)
                ! the spins of the target are the opposite of the source spins
                call intswap(spin(1), spin(2))
            else
                samplerIndex = OPP_SPIN_NO_EXCH
                ! adjust pgen
                pGen = pGen * (1.0_dp - this%pExch(ij, i_sg))
            end if
        end if
        ! get a pair of orbitals using the precomputed weights
        call this%pchb_samplers%sample(ij, samplerIndex, int(i_sg), ab, pGenHoles)
        ! split the index ab (using a table containing mapping ab -> (a,b))
        orbs = this%tgtOrbs(:, ab)
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
    end subroutine generate_double


    subroutine generate_single_GAS(nI, ilutI, nJ, ilutJ, ex_mat, tpar, store, pgen)
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel)
        integer(n_int), intent(out) :: ilutJ(0:NIfTot)
        integer, intent(out) :: ex_mat(2, maxExcit)
        logical, intent(out) :: tpar
        type(excit_gen_store_type), intent(inout), target :: store
        real(dp), intent(out) :: pGen
        @:unused_var(store)
        call gen_exc_single(general_GAS_PCHB%GASSpec, nI, ilutI, nJ, ilutJ, ex_mat, tpar, pgen)
    end subroutine

    !>  @brief
    !>  The excitation generator subroutine for PCHB.
    !>
    !>  @details
    !>  For doubles, the precomputed heat-bath weights are used.
    !>  In order to work the child classes have to override `is_allowed`
    !>  and supply a single excitation generator as function pointer in the init method.
    !>
    !>  If possible one can supply also a calc_pgen routine for single excitations.
    !>  Then it is possible to calculate the pgens for arbitrary connected configurations.
    subroutine gen_excit(this, nI, ilutI, nJ, ilutJ, ic, ex, tpar, store, pgen)
        class(GAS_PCHB_excit_gen_t), intent(in) :: this
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ic, ex(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:NIfTot)
        logical, intent(out) :: tpar
        type(excit_gen_store_type), intent(inout), target :: store
        real(dp), intent(out) :: pGen

        if (genrand_real2_dSFMT() < pSingles) then
            ic = 1
            ! defaults to uniform singles, but can be set to other excitgens
            call this%generate_single(nI, ilutI, nJ, ilutJ, ex, tpar, store, pgen)
            pgen = pgen * pSingles
        else
            ic = 2
            ! use precomputed weights to generate doubles
            call this%generate_double(nI, ilutI, nJ, ilutJ, ex, tpar, pgen)
            pgen = pgen * (1.0 - pSingles)

            if (IsNullDet(nJ)) then
                nInvalidExcits = nInvalidExcits + 1
            else
                nValidExcits = nValidExcits + 1
            end if
        end if
    end subroutine gen_excit

    !>  @brief
    !>  Deallocate the sampler and the mapping ab -> (a,b)
    subroutine finalize(this)
        class(GAS_PCHB_excit_gen_t), intent(inout) :: this
        integer :: samplerIndex
        integer(int64) :: i_sg

        call this%pchb_samplers%finalize()
        deallocate(this%tgtOrbs)
        deallocate(this%pExch)
    end subroutine

    !> @brief
    !> Placeholder function that should fail at runtime.
    real(dp) function calc_pgen_single_todo(nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2)
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(in) :: ex(2, 2), ic
        integer, intent(in) :: ClassCount2(ScratchSize), ClassCountUnocc2(ScratchSize)
        character(*), parameter :: this_routine = 'calc_pgen_single_todo'

        unused_var(nI); unused_var(ilutI); unused_var(ex); unused_var(ic);
        unused_var(ClassCount2); unused_var(ClassCountUnocc2)
        calc_pgen_single_todo = 0._dp

        call stop_all(this_routine, 'calc_pgen_single has to be implemented.')
    end function

    !>  @brief
    !>  Calculate the probability of generating a given excitation with the pchb excitgen
    !>
    !>  @param[in] nI  determinant to start from
    !>  @param[in] ex  2x2 excitation matrix
    !>  @param[in] ic  excitation level
    !>  @param[in] ClassCount2  symmetry information of the determinant
    !>  @param[in] ClassCountUnocc2  symmetry information of the virtual orbitals
    !>
    !>  @return pGen  probability of drawing this excitation with the pchb excitgen
    function calc_pgen_pchb(this, nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2) result(pgen)
        class(GAS_PCHB_excit_gen_t), intent(in) :: this
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(in) :: ex(2, 2), ic
        integer, intent(in) :: ClassCount2(ScratchSize), ClassCountUnocc2(ScratchSize)
        real(dp) :: pgen

        if (ic == 1) then
            pgen = pSingles * this%calc_pgen_single(nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2)
        else if (ic == 2) then
            pGen = (1.0 - pSingles) * this%calc_double_pgen_pchb(nI, ex)
        else
            pgen = 0.0_dp
        end if
    end function calc_pgen_pchb


    !>  @brief
    !>  Calculate the probability of drawing a given double excitation ex
    !>
    !>  @param[in] ex  2x2 excitation matrix
    !>
    !>  @return pgen  probability of generating this double with the pchb double excitgen
    function calc_double_pgen_pchb(this, nI, ex) result(pgen)
        class(GAS_PCHB_excit_gen_t), intent(in) :: this
        integer, intent(in) :: nI(nEl)
        integer, intent(in) :: ex(2, 2)
        real(dp) :: pgen
        integer :: ab, ij, nex(2, 2), samplerIndex, i_sg

        i_sg = this%indexer%idx_nI(nI)
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
                pgen = pgen * (1 - this%pExch(ij, i_sg))
            else
                ! opp spin case with exchange
                samplerIndex = OPP_SPIN_EXCH
                pgen = pgen * this%pExch(ij, i_sg)
            end if
        end if

        ! look up the probability for this excitation in the sampler
        ab = fuseIndex(nex(2, 1), nex(2, 2))
        pgen = pgen * this%pchb_samplers%get_prob(ij, samplerIndex, i_sg, ab)

    end function calc_double_pgen_pchb

end module gasci_general_pchb
