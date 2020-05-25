#include "macros.h"
module guga_pchb_excitgen

    use aliasSampling, only: aliasSamplerArray_t
    use constants, only: n_int, dp, maxExcit, int64, iout
    use bit_rep_data, only: IlutBits, GugaBits
    use SystemData, only: nel, G1, current_stepvector, t_pchb_weighted_singles, &
                          nBasis, nSpatOrbs, ElecPairs
    use FciMCData, only: excit_gen_store_type, pSingles, pDoubles
    use guga_data, only: tNewDet, ExcitationInformation_t, gen_type, excit_type
    use guga_bitrepops, only: convert_ilut_toGUGA, isProperCSF_ilut
    use dSFMT_interface, only: genrand_real2_dSFMT
    use util_mod, only: near_zero, fuseIndex, intswap
    use CalcData, only: t_matele_cutoff, matele_cutoff
    use sym_general_mod, only: ClassCountInd
    use SymExcitDataMod, only: OrbClassCount, SymLabelCounts2, &
                               sym_label_list_spat, SpinOrbSymLabel
    use UMatCache, only: gtID
    use guga_excitations, only: assign_excitinfo_values_single, &
                                createStochasticExcitation_single, &
                                pick_elec_pair_uniform_guga, &
                                excitationIdentifier_double, get_guga_integral_contrib, &
                                calc_pgen_mol_guga_single, get_excit_level_from_excitInfo
    use guga_procedure_pointers, only: gen_single_excit_guga, gen_double_excit_guga
    use guga_bitrepops, only: identify_excitation, encode_excit_info
    use bit_reps, only: decode_bit_det

    implicit none

    private

    public :: pick_orbitals_double_pchb, pick_orbitals_pure_uniform_singles, &
              calc_orbital_pgen_contr_pchb, calc_orbital_pgen_contr_start_pchb, &
              calc_orbital_pgen_contr_end_pchb, init_guga_pchb_excitgen, &
              calc_pgen_guga_pchb

    ! start with one sampler for now!
    type(aliasSamplerArray_t) :: pchb_samplers(1)
    integer, allocatable :: tgtOrbs(:,:)

    interface calc_orb_pgen_uniform_singles
        module procedure calc_orb_pgen_uniform_singles_exmat
        module procedure calc_orb_pgen_uniform_singles_excitInfo
    end interface calc_orb_pgen_uniform_singles

    interface calc_orb_pgen_guga_pchb_double
        module procedure calc_orb_pgen_guga_pchb_double_exmat
        module procedure calc_orb_pgen_guga_pchb_double_excitInfo
    end interface calc_orb_pgen_guga_pchb_double

contains

    subroutine init_guga_pchb_excitgen
        debug_function_name("init_guga_pchb_excitgen")
        integer :: ab, a, b, abMax, aerr
        integer(int64) :: memCost
        integer :: samplerIndex = 1

        write(iout,*) "Allocating GUGA PCHB excitation generator objects"
        ! total memory cost
        memCost = 0_int64
        ! initialize the mapping ab -> (a,b)
        abMax = fuseIndex(nSpatOrbs,nSpatOrbs)
        allocate(tgtOrbs(2,0:abMax), stat = aerr)
        do a = 1, nSpatOrbs
         do b = a, nSpatOrbs
            ab = fuseIndex(a,b)
            tgtOrbs(1,ab) = a
            tgtOrbs(2,ab) = b
         end do
        end do

        ! enable catching exceptions
        tgtOrbs(:,0) = 0

        call setup_pchb_sampler_conditional()

        write(iout,*) "Finished GUGA PCHB excitation generator initialization"

    contains

        subroutine setup_pchb_sampler_conditional()
            debug_function_name("setup_pchb_sampler_conditional")
            integer :: i, j, ij, ijMax
            integer(int64), allocatable :: excit_info(:)
            real(dp), allocatable :: w(:)
            ! just to be save also code up a setup function with the
            ! necessary if statements, so i can check if my optimized
            ! sampler is correct (and if it is even faster..)

            ijMax = fuseIndex(nSpatOrbs, nSpatOrbs)
            memCost = memCost + abMax*ijMax*24*2
            write(iout,*) "Excitation generator requires", &
                real(memCost,dp)/2.0_dp**30, "GB of memory"
            write(iout,*) "Generating samplers for PCHB excitation generator"

            ! weights per pair
            allocate(w(abMax), stat = aerr)
            ! excitation information per pair:
            allocate(excit_info(abMax), stat = aerr)

            ! allocate: all samplers have the same size
            call pchb_samplers(1)%setupSamplerArray(int(ijMax,int64), &
                int(abMax,int64))
            ! todo: do the same for the excit_info array!

            ! the encode_excit_info function encodes as: E_{ai}E_{bj)
            ! but for some index combinations we have to change the input so
            ! actually E_{aj}E_{bi} is encoded! convention to use only
            ! certain type of excit-types
            do i = 1, nSpatOrbs
                do j = i, nSpatOrbs
                    w = 0.0_dp
                    ij = fuseIndex(i,j)
                    excit_info = 0_int64
                    do a = 1, nSpatOrbs
                        do b = a, nSpatOrbs
                            ab = fuseIndex(a,b)
                            if (i == j) then
                                if (a == b) then
                                    ! here we only have a contribution if
                                    ! a != i
                                    if (a < i) then
                                        ! _RR_(a) -> ^RR^(i)
                                        w(ab) = get_guga_integral_contrib([i,j],a,b)
                                        excit_info(ab) = encode_excit_info(&
                                            typ = excit_type%fullstart_stop_alike, &
                                            a = a, i = i, b = b, j = j)
                                    else if (a > i) then
                                        ! _LL_(i) > ^LL^(a)
                                        w(ab) = get_guga_integral_contrib([i,j],a,b)
                                        excit_info(ab) = encode_excit_info(&
                                            typ = excit_type%fullstart_stop_alike, &
                                            a = a, i = i, b = b, j = j)
                                    end if
                                elseif (a /= b) then
                                    ! here we have to determine where (a) and
                                    ! (b) are relative to (i=j) and a == i or
                                    ! b == i are NOT allowed!
                                    if (a < i) then
                                        if (b < i) then
                                            ! _R(a) -> _RR(b) -> ^RR^(i)
                                            w(ab) = get_guga_integral_contrib([i,j],a,b)
                                            excit_info(ab) = encode_excit_info(&
                                                typ = excit_type%fullstop_raising, &
                                                a = a, i = i, b = b, j = j)
                                        else if (b > i) then
                                            ! _R(a) -> ^RL_(i) -> ^L(b)
                                            w(ab) = get_guga_integral_contrib([i,j],a,b)
                                            excit_info(ab) = encode_excit_info(&
                                                typ = excit_type%single_overlap_R_to_L, &
                                                a = a, i = i, b = b, j = j)
                                        end if
                                    else if (a > i) then
                                        ! since b > a ensured only:
                                        ! _LL_(i) -> ^LL(a) -> ^L(b)
                                        w(ab) = get_guga_integral_contrib([i,j],a,b)
                                        excit_info(ab) = encode_excit_info(&
                                            typ = excit_type%fullstart_lowering, &
                                            a = a, i = i, b = b, j = j)
                                    end if
                                end if
                            else if ( i /= j) then
                                if (a == b) then
                                    ! a == i or a == j NOT allowed!
                                    if (a < i) then
                                        ! _RR_(a) -> ^RR(i) -> ^R(j)
                                        w(ab) = get_guga_integral_contrib([i,j],a,b)
                                        excit_info(ab) = encode_excit_info(&
                                            typ = excit_type%fullstart_raising, &
                                            a = a, i = i, b = b, j = j)
                                    else if (a > i .and. a < j) then
                                        ! _L(i) -> ^LR_(a) -> ^R(j)
                                        w(ab) = get_guga_integral_contrib([i,j],a,b)
                                        excit_info(ab) = encode_excit_info(&
                                            typ = excit_type%single_overlap_L_to_R, &
                                            a = a, i = i, b = b, j = j)
                                    else if (a > j) then
                                        ! _L(i) -> _LL(j) -> ^LL^(a)
                                        w(ab) = get_guga_integral_contrib([i,j],a,b)
                                        excit_info(ab) = encode_excit_info(&
                                            typ = excit_type%fullstop_lowering, &
                                            a = a, i = i, b = b, j = j)
                                    end if
                                else if (a /= b) then
                                    ! this is the most general case. a lot of IFs
                                    if (a < i) then
                                        if (b < i) then
                                            w(ab) = get_guga_integral_contrib([i,j],a,b)
                                            ! in E_{ai}E_{bi} form this would be
                                            ! _R(a) -> RR_(b) -> ^RR(i) -> ^R(j)
                                            ! which would have the opposite
                                            ! sign conventions for the x1
                                            ! elements as in the Shavitt 81
                                            ! paper.
                                            ! (technically it would not matter
                                            !  here since both are flipped,
                                            !  but lets stay consistent!)
                                            ! for E_{bj}E_{ai}
                                            ! _R(a) -> _RR(b) -> RR^(i) -> ^R(j)
                                            ! it would fit! so:
                                            excit_info(ab) = encode_excit_info(&
                                                typ = excit_type%double_raising, &
                                                a = b, i = j, b = a, j = i)
                                        ! else if (b == i) then
                                            ! b == i also NOT allowed here,
                                            ! since this would correspond to
                                            ! a single!
                                        else if (b > i .and. b < j) then
                                            w(ab) = get_guga_integral_contrib([i,j],a,b)
                                            ! for E_{ai}E_{bj} this would
                                            ! correspond to a non-overlap:
                                            ! _R(a) -> ^R(i) + _R(b) > ^R(j)
                                            ! which are not directly sampled
                                            ! However for E_{aj}E_{bj} this is
                                            ! _R{a} -> _LR(i) -> ^LR(b) -> ^R(j)
                                            excit_info(ab) = encode_excit_info(&
                                                typ = excit_type%double_R_to_L_to_R, &
                                                a = a, i = j, b = b, j = i)
                                        else if (b == j) then
                                            w(ab) = get_guga_integral_contrib([i,j],a,b)
                                            ! here we also have to switch
                                            ! indices to E_{aj}E_{bi} to get:
                                            ! _R(i) -> _LR(a) -> ^LR^(j)
                                            excit_info(ab) = encode_excit_info(&
                                                typ = excit_type%fullstop_R_to_L, &
                                                a = a, i = j, b = b, j = i)

                                        else if (b > j) then
                                            w(ab) = get_guga_integral_contrib([i,j],a,b)
                                            ! here we have to switch to
                                            ! E_{aj}E_{bi} to get:
                                            ! _R(a) > _LR(i) -> LR^(j) -> ^L(b)
                                            excit_info(ab) = encode_excit_info(&
                                                typ = excit_type%double_R_to_L, &
                                                a = a, i = j, b = b, j = i)
                                        end if
                                    else if (a == i) then
                                        ! b > i is ensured here since b > a in here!
                                        if (b < j) then
                                            w(ab) = get_guga_integral_contrib([i,j],a,b)
                                            ! here we have to switch again:
                                            ! E_{aj}E_{bi}:
                                            ! _RL_(i) -> ^LR(b) -> ^R(j)
                                            excit_info(ab) = encode_excit_info(&
                                                typ = excit_type%fullstart_L_to_R, &
                                                a = a, i = j, b = b, j = i)
                                        else if (b == j) then
                                            w(ab) = get_guga_integral_contrib([i,j],a,b)
                                            ! switch again: E_{aj}E_{bi}
                                            ! _RL_(i) -> _RL_(j)
                                            excit_info(ab) = encode_excit_info(&
                                                typ = excit_type%fullstart_stop_mixed, &
                                                a = a, i = j, b = b, j = i)
                                        else if (b > j) then
                                            w(ab) = get_guga_integral_contrib([i,j],a,b)
                                            ! switch again: E_{aj}E_{bi}
                                            ! _RL_(i) -> ^RL(j) -> ^L(b)
                                            excit_info(ab) = encode_excit_info(&
                                                typ = excit_type%fullstart_R_to_L, &
                                                a = a, i = j, b = b, j = i)
                                        end if
                                    else if (a > i .and. a < j) then
                                        ! b > a still ensured!
                                        if (b < j) then
                                            w(ab) = get_guga_integral_contrib([i,j],a,b)
                                            ! switch: E_{aj}E_{bi}:
                                            ! _L(j) -> _RL(a) -> ^LR(b) -> ^R(j)
                                            excit_info(ab) = encode_excit_info(&
                                                typ = excit_type%double_L_to_R, &
                                                a = a, i = j, b = b, j = i)

                                        else if (b == j) then
                                            w(ab) = get_guga_integral_contrib([i,j],a,b)
                                            ! switch: E_{aj}E_{bi}
                                            ! _L(i) -> _RL(a) -> ^RL^(j)
                                            excit_info(ab) = encode_excit_info(&
                                                typ = excit_type%fullstop_L_to_R, &
                                                a = a, i = j, b = b, j = i)
                                        else if (b > j) then
                                            w(ab) = get_guga_integral_contrib([i,j],a,b)
                                            ! switch: E_{aj}E_{bi}
                                            ! _L(i) -> _RL(a) - > ^RL(j) -> ^L(b)
                                            excit_info(ab) = encode_excit_info(&
                                                typ = excit_type%double_L_to_R_to_L, &
                                                a = a, i = j, b = b, j = i)
                                        end if
                                    ! else if (a == j) then
                                        ! a == j also NOT allowed here!
                                    else if (a > j) then
                                        ! b > a > j implied here!
                                        w(ab) = get_guga_integral_contrib([i,j],a,b)
                                        ! E_{ai}E_{bj} would lead to:
                                        ! _L(i) -> LL_(j) -> ^LL(a) -> ^L(b)
                                        ! which has the correct sign convention
                                        ! as in the Shavitt 81 paper
                                        excit_info(ab) = encode_excit_info(&
                                            typ = excit_type%double_lowering, &
                                            a = a, i = i, b = b, j = j)
                                    end if
                                end if
                            end if
                        end do
                    end do

                    call pchb_samplers(samplerIndex)%setupEntry(ij,w)
                    ! todo: do the same for the excit_info array!
                end do
            end do



        end subroutine setup_pchb_sampler_conditional

!         subroutine setup_optimized_pchb_sampler()
!             debug_function_name("setup_optimized_pchb_sampler")
!             integer :: i, j, ij, ijMax
!             integer(int64), allocatable :: excit_info(:)
!             real(dp), allocatable :: w(:)
!
!             ijMax = fuseIndex(nSpatOrbs, nSpatOrbs)
!             memCost = memCost + abMax*ijMax*24*2
!             write(iout,*) "Excitation generator requires", &
!                 real(memCost,dp)/2.0_dp**30, "GB of memory"
!             write(iout,*) "Generating samplers for PCHB excitation generator"
!
!             ! weights per pair
!             allocate(w(abMax), stat = aerr)
!             ! excitation information per pair:
!             allocate(excit_info(abMax), stat = aerr)
!
!             ! allocate: all samplers have the same size
!             call pchb_samplers(1)%setupSamplerArray(int(ijMax,int64), &
!                 int(abMax,int64))
!             ! todo: do the same for the excit_info array!
!
!             ! do the loops in an optimized way to avoid IF statements
!             do i = 1, nSpatOrbs
!                 ! deal with i = j:
!                 w = 0.0_dp
!                 j = i
!                 ! a -> i
!                 do a = 1, i - 1 ! a = i = j not valid
!                     ! deal with b = a
!                     b = a
!
!                     ! we know here we have a RR(a) -> RR(i)
!
!                     ab = fuseIndex(a,b)
!                     w(ab) = get_guga_integral_contrib([i,j], a, b)
!                     ! the excitation is encoded like: e_{ij,kl}
!                     excit_info(ab) = encode_excit_info( &
!                         typ = excit_type%fullstart_stop_alike, &
!                         i = a, j = i, k = b, l = j)
!
!
!                     ! now deal with the rest of b:
!                     do b = a + 1, i - 1 ! b = i = j also not valid!
!
!                         ! here we have _R(a) -> _RR(b) -> ^RR^(ij)
!
!                         ab = fuseIndex(a,b)
!                         w(ab) = get_guga_integral_contrib([i,j], a, b)
!                         excit_info(ab) = encode_excit_info( &
!                             typ = excit_type%fullstop_raising, &
!                             i = a, j = i, k = b, l = j)
!
!                     end do
!
!                     do b = i + 1, nSpatOrbs
!                         ! here we have: single overlap raising into lowering
!                         ab = fuseIndex(a,b)
!                         w(ab) = get_guga_integral_contrib([i,j],a,b)
!                         excit_info(ab) = encode_excit_info(&
!                             typ = excit_type%single_overlap_R_to_L, &
!                             i = a, j = i, k = b, l = j)
!
!                     end do
!                 end do
!
!                 ! a -> nSpatOrbs
!                 do a = i + 1, nSpatOrbs
!                     ! again deal with b = a
!                     b = a
!                     ! todo.. the rest.. gonna be tough, but then i have so
!                     ! much information already!
!
!                 ij = fuseIndex(i,j)
!                 call pchb_samplers(samplerIndex)%setupEntry(ij,w)
!                 ! todo: do the same for the excit_info list!
!
!
!                 do
!
!
!         end subroutine setup_optimized_pchb_sampler

        subroutine setup_guga_pchb_sampler()
            debug_function_name("setup_guga_pchb_sampler")
            integer :: i, j, ij, ijMax
            real(dp), allocatable :: w(:)
            ! number of possible source orbital pairs
            ijMax = fuseIndex(nSpatOrbs,nSpatOrbs)

            memCost = memCost + abMax*ijMax*24
            write(iout,*) &
                "Excitation generator requires", real(memCost,dp)/2.0_dp**30, "GB of memory"
            write(iout,*) "Generating samplers for PCHB excitation generator"

            ! for now do it in the most simply way and do not include any
            ! GUGA restrictions here. Just check if it works and then
            ! improve!

            ! weights per pair
            allocate(w(abMax), stat = aerr)

            ! allocate: all samplers have the same size
            call pchb_samplers(samplerIndex)%setupSamplerArray(int(ijMax,int64),int(abMax,int64))

            do i = 1, nSpatOrbs
                do j = i, nSpatOrbs
                    w = 0.0_dp
                    do a = 1, nSpatOrbs
                        do b = a, nSpatOrbs
                            ab = fuseIndex(a,b)
                            w(ab) = get_guga_integral_contrib([i,j],a,b)
                        end do
                    end do
                    ij = fuseIndex(i,j)
                    call pchb_samplers(samplerIndex)%setupEntry(ij,w)
                end do
            end do

            deallocate(w)

        end subroutine setup_guga_pchb_sampler

    end subroutine init_guga_pchb_excitgen

    subroutine finalize_pchb_excitgen_guga()

        call pchb_samplers(1)%samplerArrayDestructor()
        deallocate(tgtOrbs)

    end subroutine finalize_pchb_excitgen_guga

    subroutine pick_orbitals_double_pchb(ilut, nI, excitInfo, pgen)
        debug_function_name("pick_orbitals_double_pchb")
        integer(n_int), intent(in) :: ilut(0:IlutBits%len_tot)
        integer, intent(in) :: nI(nel)
        type(ExcitationInformation_t), intent(out) :: excitInfo
        real(dp), intent(out) :: pgen

        integer :: src(2), sym_prod, sum_ml, i, j, a, b, orbs(2), ij, ab
        real(dp) :: pgen_elec, pgen_orbs

        ! maybe I will also produce a weighted electron pickin in the
        ! GUGA formalism.. but for now pick them uniformly:
        call pick_elec_pair_uniform_guga(nI, src, sym_prod, sum_ml, pgen_elec)
        ASSERT( src(1) < src(2) )

        i = gtID(src(1))
        j = gtID(src(2))

        ! use the sampler for this electron pair -> order of src electrons
        ! does not matter
        ij = fuseIndex(i, j)

        ! get a pair of orbitals using the precomputed weights
        call pchb_samplers(1)%aSample(ij,ab,pgen_orbs)

        ! unfortunately, there is a super-rare case when, due to floating point error,
        ! an excitation with pGen=0 is created. Those are invalid, too
        if(near_zero(pgen_orbs)) then
            excitInfo%valid = .false.
            ! Yes, print. Those events are signficant enough to be always noted in the output
            print *, "WARNING: Generated excitation with probability of 0"
            return
        endif

        ! split the index ab (using a table containing mapping ab -> (a,b))
        orbs = tgtOrbs(:,ab)

        a = orbs(1)
        b = orbs(2)

        ! check if the picked orbs are a valid choice - if they are the same, match one
        ! occupied orbital or are zero (maybe because there are no allowed picks for
        ! the given source) abort
        ! these are the 'easy' checks for GUGA.. more checks need to be done
        ! to see if it is actually a valid combination..
        if (any(orbs == 0) .or. (current_stepvector(a) == 3) .or. &
            (current_stepvector(b) == 3) .or. (a == b .and. &
            current_stepvector(b) /= 0)) then
            excitInfo%valid = .false.
            return
        end if

        excitInfo = excitationIdentifier_double(a, i, b, j)

        pGen = pgen_elec * pgen_orbs

    end subroutine pick_orbitals_double_pchb

    subroutine pick_orbitals_pure_uniform_singles(ilut, nI, excitInfo, pgen)
        debug_function_name("pick_orbitals_pure_uniform_singles")
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        integer, intent(in) :: nI(nel)
        type(ExcitationInformation_t), intent(out) :: excitInfo
        real(dp), intent(out) :: pgen

        integer :: elec

        unused_var(ilut)

        ASSERT(isProperCSF_ilut(ilut, .true.))

        ! init for safety and so we can return on abort
        pgen = 0.0_dp

        ! pick random electron
        ! have to modify pgen for doubly occupied orbs! since double the
        ! chance!
        elec = 1 + floor(genrand_real2_dSFMT() * nel)

        call pick_uniform_spatial_hole(nI, elec, excitInfo, pgen)

    end subroutine pick_orbitals_pure_uniform_singles

    subroutine pick_uniform_spatial_hole(nI, elec, excitInfo, pgen)
        debug_function_name("pick_uniform_spatial_hole")
        integer, intent(in) :: nI(nel), elec
        type(ExcitationInformation_t), intent(out) :: excitInfo
        real(dp), intent(out) :: pgen

        integer :: so_elec, s_elec, cc_i, nOrb, sym_index, attempts, orb, &
                   s_orb
        integer, parameter :: limit = 250

        pgen = 0.0_dp

        ! electron spin-orbital
        so_elec = nI(elec)

        ! electron spatial index
        s_elec = gtID(so_elec)

        ! get the symmetry index for this electron
        cc_i = ClassCountInd(1, SpinOrbSymLabel(so_elec), G1(so_elec)%ml)

        ! and the number of orbitals
        nOrb = OrbClassCount(cc_i)

        if (nOrb == 0) return

        ! get the symmetry index for later use
        sym_index = SymLabelCounts2(1, cc_i)

        ! now keep drawing from the symmetry orbitals until we pick an
        ! 'empty' (not doubly occupied!) one
        attempts = 0
        do while (attempts < limit)

            orb = 1 + floor(genrand_real2_dSFMT() * nOrb)

            s_orb = sym_label_list_spat(sym_index + orb - 1)

            ! if the spatial orbital is not doubly occupied an GUGA
            ! excitation (without taking into accound any other restrictions)
            ! can be possible
            ! and it must not be the original electron spatial index!
            if (current_stepvector(s_orb) /= 3 .and. s_orb /= s_elec) exit

            attempts = attempts + 1
#ifdef DEBUG_
            if (attempts > 200) then
                print *, "closing to 250 limit in random single orbital picking"
            end if
#endif
        end do

        if (s_orb < s_elec) then
            excitInfo = assign_excitinfo_values_single(gen_type%R, s_orb, s_elec, &
                s_orb, s_elec)
        else
            excitInfo = assign_excitinfo_values_single(gen_type%L, s_orb, s_elec, &
                s_elec, s_orb)
        end if

        ! do i need nOrb now or the actual number of unoccupied in the
        ! CSF? i think the second..
        pgen = 1.0_dp / real(nOrb * nel, dp)

        if (current_stepvector(s_elec) == 3) pgen = 2.0_dp * pgen

    end subroutine pick_uniform_spatial_hole

    function calc_orb_pgen_uniform_singles_exmat(ex) result(pgen)
        integer, intent(in) :: ex(2,2)
        real(dp) :: pgen
        debug_function_name("calc_orb_pgen_uniform_singles_exmat")

        integer :: src(2), so_elec, cc_i, nOrb
#ifdef DEBUG_
        integer :: tgt(2)

        tgt = get_tgt(ex)
        ASSERT(current_stepvector(gtID(tgt(1))) /= 3)
#endif
        src = get_src(ex)

        ASSERT(all(ex(:,2) == 0))
        ASSERT(all(ex(:,1) > 0))
        ASSERT(all(ex(:,1) <= nBasis))
        ASSERT(ex(1,1) /= ex(2,1))
        ASSERT(current_stepvector(gtID(src(1))) /= 0)

        so_elec = src(1)

        cc_i = ClassCountInd(1, SpinOrbSymLabel(so_elec), G1(so_elec)%ml)

        nOrb = OrbClassCount(cc_i)

        pgen = 1.0_dp / real(nOrb * nel, dp)

    end function calc_orb_pgen_uniform_singles_exmat

    function calc_orb_pgen_uniform_singles_excitInfo(excitInfo) result(pgen)
        type(ExcitationInformation_t), intent(in) :: excitInfo
        real(dp) :: pgen

        pgen = 0.0_dp

    end function calc_orb_pgen_uniform_singles_excitInfo

    function calc_pgen_guga_pchb(ilutI, ilutJ, excitInfo_in) result(pgen)
        debug_function_name("calc_pgen_guga_pchb")
        integer(n_int), intent(in) :: ilutI(0:GugaBits%len_tot), ilutJ(GugaBits%len_tot)
        type(ExcitationInformation_t), intent(in), optional :: excitInfo_in
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: pgen
        integer :: ic
        integer :: nI(nel), nJ(nel)

        if (present(excitInfo_in)) then
            excitInfo = excitInfo_in
        else
            excitInfo = identify_excitation(ilutI, ilutJ)
        end if

        ic = get_excit_level_from_excitInfo(excitInfo)

        if (ic == 1) then
            if (t_pchb_weighted_singles) then
                call decode_bit_det(nI, ilutI)
                call decode_bit_det(nJ, ilutJ)
                pgen = calc_pgen_mol_guga_single(ilutI, nI, ilutJ, &
                    nJ, excitInfo)
            else
                pgen = calc_orb_pgen_uniform_singles(excitInfo)
            end if

            pgen = pgen * pSingles

        else if (ic == 2) then

            pgen = pDoubles * calc_orb_pgen_guga_pchb_double(excitInfo)

        else
            pgen = 0.0_dp
        end if

    end function calc_pgen_guga_pchb

    function calc_orb_pgen_guga_pchb_double_exmat(ex) result(pgen)
        debug_function_name("calc_orb_pgen_guga_pchb_double_exmat")
        integer, intent(in) :: ex(2,2)
        real(dp) :: pgen
        integer :: ij, ab, nex(2,2)
        real(dp) :: p_elec

        nex = gtID(ex)
        ij = fuseIndex(nex(1,1), nex(1,2))

        p_elec = 1.0_dp / real(ElecPairs, dp)

        ab = fuseIndex(nex(2,1), nex(2,2))
        pgen = p_elec * pchb_samplers(1)%aGetProb(ij, ab)

    end function calc_orb_pgen_guga_pchb_double_exmat

    function calc_orb_pgen_guga_pchb_double_excitInfo(excitInfo) result(pgen)
        debug_function_name("calc_orb_pgen_guga_pchb_double_excitInfo")
        type(ExcitationInformation_t), intent(in) :: excitInfo
        real(dp) :: pgen

        pgen = 0.0_dp


    end function calc_orb_pgen_guga_pchb_double_excitInfo

    ! I need the pgen-recalculation routines for exchange type excitations
    ! also for the PCHB excit-gen
    subroutine calc_orbital_pgen_contr_pchb(ilut, occ_orbs, cpt_a, cpt_b)
        debug_function_name("calc_orbital_pgen_contr_pchb")
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        integer, intent(in) :: occ_orbs(2)
        real(dp), intent(out) :: cpt_a, cpt_b

        integer :: i, j, ij
        unused_var(ilut)

        ! this function can in theory be called with both i < j and i > j..
        ! to take the correct values here!
        i = minval(gtID(occ_orbs))
        j = maxval(gtID(occ_orbs))

        ! i think i have to consider both the above and below contribution..
        ! but i am not so sure how.. in this PCHB case..
        ! maybe in PCHB those 2 are just the same.. i am confused
        ij = fuseIndex(i, j)

        ! and both I and J are electron and hole indices here
        cpt_a = pchb_samplers(1)%aGetProb(ij, ij)

        ! and in the PCHB there is no difference between the 2!
        cpt_b = cpt_a

    end subroutine calc_orbital_pgen_contr_pchb


    ! i think it would be better if i 'just' reimplement:
    function calc_orbital_pgen_contr_start_pchb(occ_orbs, a) result(orb_pgen)
        debug_function_name("calc_orbital_pgen_contr_start_pchb")
        integer, intent(in) :: occ_orbs(2), a
        real(dp) :: orb_pgen

        integer :: i, j, ij, ab

        ! i and j should here always be provided with i < j
        ASSERT( i < j )

        ! depending on type (R->L / L->R) a can be > j or < j, but always > i
        ASSERT( i < a )
        i = gtID(occ_orbs(1))
        j = gtID(occ_orbs(2))

        ! here i is both electron and hole index!

        ij = fuseIndex(i, j)
        ab = fuseIndex(i, a)

        orb_pgen = pchb_samplers(1)%aGetProb(ij, ab)

    end function calc_orbital_pgen_contr_start_pchb

    function calc_orbital_pgen_contr_end_pchb(occ_orbs, a) result(orb_pgen)
        debug_function_name("calc_orbital_pgen_contr_end_pchb")
        integer, intent(in) :: occ_orbs(2), a
        real(dp) :: orb_pgen

        integer :: i, j, ij, ab

        ! I can assert i < j here.. as this is how this function should
        ! be called!
        ASSERT( i < j )
        i = gtID(occ_orbs(1))
        j = gtID(occ_orbs(2))

        ! and j is at the same time electron and hole index!

        ! depending on L->R or R->L type a can be > i o r < i, but always < j!
        ASSERT( a < j )

        ! j here is both elec and hole ind!
        ij = fuseIndex(i,j)

        ab = fuseIndex(a,j)

        orb_pgen = pchb_samplers(1)%aGetProb(ij,ab)

    end function calc_orbital_pgen_contr_end_pchb

end module guga_pchb_excitgen

