#include "macros.h"

module excit_gen_5

    use excit_gens_int_weighted, only: gen_single_4ind_ex, pgen_single_4ind, &
                                       get_paired_cc_ind, select_orb, &
                                       opp_spin_pair_contrib, &
                                       same_spin_pair_contrib, &
                                       pick_biased_elecs, &
                                       pick_weighted_elecs, select_orb, &
                                       pgen_select_orb, pgen_weighted_elecs
    use SymExcitDataMod, only: SpinOrbSymLabel, SymInvLabel, ScratchSize
    use FciMCData, only: excit_gen_store_type, pSingles, pDoubles
    use SystemData, only: G1, tUHF, tStoreSpinOrbs, nbasis, nel, &
                          tGen_4ind_part_exact
    use SymExcit3, only: CountExcitations3, GenExcitations3
    use GenRandSymExcitNUMod, only: init_excit_gen_store
    use DetBitOps, only: ilut_lt, ilut_gt, EncodeBitDet
    use Determinants, only: write_det, get_helement
    use dSFMT_interface, only: genrand_real2_dSFMT
    use procedure_pointers, only: get_umat_el
    use sym_general_mod, only: ClassCountInd
    use bit_rep_data, only: NIfTot, NIfD
    use bit_reps, only: decode_bit_det
    use get_excit, only: make_double
    use UMatCache, only: gtid
    use constants
    use sort_mod
    use util_mod
    implicit none

contains

    subroutine gen_excit_4ind_weighted2 (nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                     ExcitMat, tParity, pGen, HelGen, store)

        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), IC, ExcitMat(2,2)
        logical, intent(out) :: tParity
        real(dp), intent(out) :: pGen
        HElement_t, intent(out) :: HElGen
        type(excit_gen_store_type), intent(inout), target :: store
        integer(n_int), intent(out) :: ilutJ(0:NIfTot)
        character(*), parameter :: this_routine = 'gen_excit_4ind_weighted'

        real(dp) :: pgen2

        HElGen = HEl_zero

        ! Choose if we want to do a single or a double excitation
        ! TODO: We can (in principle) re-use this random number by subdivision
        if (genrand_real2_dSFMT() < pSingles) then

            ic = 1
            call gen_single_4ind_ex (nI, ilutI, nJ, ilutJ, ExcitMat, &
                                     tParity, pGen)
            pgen = pgen * pSingles

        else

            ! OK, we want to do a double excitation
            ic = 2
            call gen_double_4ind_ex2 (nI, ilutI, nJ, ilutJ, ExcitMat, tParity, &
                                      pGen)
            pgen = pgen * pDoubles

        end if

        ! And a careful check!
#ifdef __DEBUG
        if (.not. IsNullDet(nJ)) then
            pgen2 = calc_pgen_4ind_weighted2(nI, ilutI, ExcitMat, ic)
            if (abs(pgen - pgen2) > 1.0e-6_dp) then
                write(6,*) 'Calculated and actual pgens differ.'
                write(6,*) 'This will break HPHF calculations'
                call write_det(6, nI, .false.)
                write(6, '(" --> ")', advance='no')
                call write_det(6, nJ, .true.)
                write(6,*) 'Excitation matrix: ', ExcitMat(1,1:ic), '-->', &
                           ExcitMat(2,1:ic)
                write(6,*) 'Generated pGen:  ', pgen
                write(6,*) 'Calculated pGen: ', pgen2
                call stop_all(this_routine, "Invalid pGen")
            end if
        end if
#endif

    end subroutine


    function calc_pgen_4ind_weighted2 (nI, ilutI, ex, ic) &
            result(pgen)

        ! What is the probability of the excitation _from_ determinant nI
        ! described by the excitation matrix ex, and the excitation level ic,
        ! being generated according to the 4ind_weighted excitaiton generator?

        integer, intent(in) :: nI(nel), ex(2,2), ic
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        real(dp) :: pgen
        character(*), parameter :: this_routine = 'calc_pgen_4ind_weighted'

        integer :: iSpn, src(2), tgt(2)
        real(dp) :: cum_sums(2), int_cpt(2), cpt_pair(2), sum_pair(2)

        if (ic == 1) then

            ! Singles are calculated in the same way for the _ex and the
            ! _reverse excitation generators
            pgen = pSingles * pgen_single_4ind (nI, ilutI, ex(1,1), ex(2,1))

        else if (ic == 2) then

            ! This is a double excitation...
            pgen = pDoubles

            ! Get properties of source electrons
            src = ex(1, :)
            if (is_beta(src(1)) .eqv. is_beta(src(2))) then
                if (is_beta(src(1))) then
                    iSpn = 1
                else
                    iSpn = 3
                end if
            else
                iSpn = 2
            end if

            ! Select a pair of electrons in a weighted fashion
            pgen = pgen * pgen_weighted_elecs(nI, src)

            ! Obtain the probability components of picking the electrons in
            ! either A--B or B--A order
            tgt = ex(2, :)
            call pgen_select_a_orb(ilutI, src, tgt(1), iSpn, int_cpt(1), &
                                   cum_sums(1))
            call pgen_select_orb(ilutI, src, tgt(1), tgt(2), int_cpt(2), &
                                 cum_sums(2))

            call pgen_select_a_orb(ilutI, src, tgt(2), iSpn, cpt_pair(1), &
                                   sum_pair(1))
            call pgen_select_orb(ilutI, src, tgt(2), tgt(1), cpt_pair(2), &
                                 sum_pair(2))

            ! And adjust the probability for the components
            pgen = pgen * (product(int_cpt) / product(cum_sums) + &
                           product(cpt_pair) / product(sum_pair))

        else

            ! Deal with some outsider cases that can leak through the HPHF
            ! system.
            pgen = 0.0_dp

        end if

    end function


    subroutine gen_double_4ind_ex2 (nI, ilutI, nJ, ilutJ, ex, par, pgen)

        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ex(2,2)
        integer(n_int), intent(out) :: ilutJ(0:NIfTot)
        logical, intent(out) :: par
        real(dp), intent(out) :: pgen
        character(*), parameter :: this_routine = 'gen_double_4ind_ex2'

        integer :: elecs(2), src(2), orbs(2), ispn, sum_ml, cc_a, cc_b
        integer :: sym_product
        real(dp) :: int_cpt(2), cum_sum(2), cpt_pair(2), sum_pair(2)

        ! Pick the electrons in a weighted fashion
        call pick_weighted_elecs(nI, elecs, src, sym_product, ispn, sum_ml, &
                                 pgen)
        !call pick_biased_elecs(nI, elecs, src, sym_product, ispn, sum_ml, pgen)

        ! Select the A orbital _excluding_ knowledge of symmetry information.
        ! Only exclude terms that have no coupling elements.
        ! TODO: The pick_a_orb selection should only work if there are avail
        !       b orbs of the appropriate symmetry!
        orbs(1) = pick_a_orb(ilutI, src, iSpn, int_cpt(1), cum_sum(1))

        ! Select the B orbital, in the same way as before!!
        ! The symmetry of this second orbital depends on that of the first.
        if (orbs(1) /= 0) then
            cc_a = ClasSCountInd(orbs(1))
            cc_b = get_paired_cc_ind(cc_a, sym_product, sum_ml, iSpn)
            orbs(2) = select_orb (ilutI, src, cc_b, orbs(1), int_cpt(2), &
                                  cum_sum(2))
        end if

        if (any(orbs == 0)) then
            nJ(1) = 0
            return
        end if

        ! Calculate the pgens. Note that all of these excitations can be
        ! selected as both A--B or B--A. So these need to be calculated
        ! explicitly.
        ASSERT(tGen_4ind_part_exact)
        call pgen_select_a_orb(ilutI, src, orbs(2), iSpn, cpt_pair(1), &
                               sum_pair(1))
        call pgen_select_orb(ilutI, src, orbs(2), orbs(1), &
                             cpt_pair(2), sum_pair(2))
        pgen = pgen * (product(int_cpt) / product(cum_sum) + &
                       product(cpt_pair) / product(sum_pair))

        ! And generate the actual excitation
        call make_double (nI, nJ, elecs(1), elecs(2), orbs(1), orbs(2), &
                          ex, par)
        ilutJ = ilutI
        clr_orb (ilutJ, src(1))
        clr_orb (ilutJ, src(2))
        set_orb (ilutJ, orbs(1))
        set_orb (ilutJ, orbs(2))

    end subroutine gen_double_4ind_ex2


    function pick_a_orb(ilut, src, ispn, cpt, cum_sum) result(orb)

        integer(n_int), intent(in) :: ilut(0:NifTot)
        integer, intent(in) :: src(2), iSpn
        real(dp), intent(out) :: cpt, cum_sum
        integer :: orb
        character(*), parameter :: t_r = 'pick_a_orb'
        character(*), parameter :: this_routine = t_r

        real(dp) :: cum_arr(nbasis), r, cum_tst, cpt_tst
        integer :: start_ind, srcid(2)
        logical :: occa, occb

        logical :: beta, parallel, valid

        ! Just in case. eeep.
        ! Note that we scale spin/spatial orbitals here with factors of 2
        ! which need more care if using spin orbitals
        if (tUHF .or. tStoreSpinOrbs) &
            call stop_all(this_routine, "ASSUMES RHF orbitals")

        if (iSpn == 1) then
            parallel = .true.
            beta = .true.
        else if (iSpn == 2) then
            parallel = .false.
            beta = .true.
        else
            parallel = .true.
            beta = .false.
        end if

        cum_sum = 0
        do orb = 1, nbasis

            ! TODO: This should be doable in a better way...
            if (beta .eqv. is_beta(src(1))) then
                srcid = gtID(src)
            else
                ASSERT(iSpn == 2)
                ASSERT(beta .eqv. is_beta(src(2)))
                srcid(1) = gtID(src(2))
                srcid(2) = gtID(src(1))
            end if

            valid = .true.
            !if (.not. parallel .or. is_beta(orb) .neqv. beta) valid = .false.

            if (IsOcc(ilut, orb)) valid = .false.
            if (parallel .and. (is_beta(orb) .neqv. beta)) valid = .false.

            cpt = 0
            if (valid) then
                ! Get the correct element depending on spin terms
                if (parallel) then
                    cpt = same_spin_pair_contrib(srcid(1), srcid(2), orb, -1)
                else
                    cpt = opp_spin_pair_contrib(srcid(1), srcid(2), orb, -1)
                end if
            end if

            cum_sum = cum_sum + cpt
            cum_arr(orb) = cum_sum
        end do

        ! And exit if invalid
        if (cum_sum == 0) then
            orb = 0
            return
        end if

        ! Pick the orbital, and extract the relevant list components for
        ! probability generation purposes
        r = genrand_real2_dSFMT() * cum_sum
        orb = binary_search_first_ge(cum_arr, r)
        if (orb == 1) then
            cpt = cum_arr(1)
        else
            cpt = cum_arr(orb) - cum_arr(orb - 1)
        end if

#ifdef __DEBUG
        call pgen_select_a_orb(ilut, src, orb, iSpn, cpt_tst, cum_tst)
        if (abs(cpt_tst - cpt) > 1e-6 .or. abs(cum_tst - cum_sum) > 1e-6) then
            call stop_all(t_r, 'Calculated probability does not match')
        end if
#endif

    end function pick_a_orb

    subroutine pgen_select_a_orb(ilut, src, orb_in, iSpn, cpt_out, cum_sum)

        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer, intent(in) :: src(2), orb_in, iSpn
        real(dp), intent(out) :: cpt_out, cum_sum
        character(*), parameter :: t_r = 'pgen_select_a_orb'
        character(*), parameter :: this_routine = t_r

        logical :: beta, parallel, valid
        integer :: srcid(2), orb
        real(dp) :: cpt

        ASSERT(is_beta(src(1)) .or. iSpn /= 1)
        ASSERT(is_beta(src(2)) .or. iSpn /= 1)
        ASSERT(.not. is_beta(src(1)) .or. iSpn /= 3)
        ASSERT(.not. is_beta(src(2)) .or. iSpn /= 3)

        if (iSpn == 1) then
            parallel = .true.
            beta = .true.
        else if (iSpn == 2) then
            parallel = .false.
            beta = .true.
        else
            parallel = .true.
            beta = .false.
        end if

        cum_sum = 0
        do orb = 1, nbasis

            ! TODO: This should be doable in a better way...
            if (beta .eqv. is_beta(src(1))) then
                srcid = gtID(src)
            else
                ASSERT(iSpn == 2)
                srcid(1) = gtID(src(2))
                srcid(2) = gtID(src(1))
            end if

            valid = .true.
            if (IsOcc(ilut, orb)) valid = .false.
            if (parallel .and. (is_beta(orb) .neqv. beta)) valid = .false.

            cpt = 0
            if (valid) then
                if (parallel) then
                    cpt = same_spin_pair_contrib(srcid(1), srcid(2), orb, -1)
                else
                    cpt = opp_spin_pair_contrib(srcid(1), srcid(2), orb, -1)
                end if
            end if

            if (orb_in == orb) cpt_out = cpt
            cum_sum = cum_sum + cpt
        end do

        ! For us to be calculating the likelihood of selecting the B electron,
        ! there must have been a possibility of selecting the A electron
        ! originally.
        if (cum_sum == 0) then
            call stop_all(t_r, 'Invalid cumulative sum')
        end if

    end subroutine

    subroutine test_excit_gen_take2 (ilut, iterations)

        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer, intent(in) :: iterations
        character(*), parameter :: this_routine = 'test_excit_gen_take2'

        integer :: src_det(nel), det(nel), nsing, ndoub, nexcit, ndet, ex(2,2)
        integer :: flag, ngen, pos, iunit, i, ic
        type(excit_gen_store_type) :: store
        integer(n_int) :: tgt_ilut(0:NifTot)
        integer(n_int), allocatable :: det_list(:,:)
        real(dp), allocatable :: contrib_list(:)
        logical, allocatable :: generated_list(:)
        logical :: found_all, par
        real(dp) :: contrib, pgen
        HElement_t :: helgen, hel

        ! Decode the determiant
        call decode_bit_det (src_det, ilut)

        ! Initialise
        call init_excit_gen_store (store)

        ! How many connected determinants are we expecting?
        call CountExcitations3 (src_det, 2, nsing, ndoub)
        nexcit = nsing + ndoub
        allocate(det_list(0:NIfTot, nexcit))

        ! Loop through all of the possible excitations
        ndet = 0
        found_all = .false.
        ex = 0
        flag = 2
        write(6,'("*****************************************")')
        write(6,'("Enumerating excitations")')
        write(6,'("Starting from: ")', advance='no')
        call write_det (6, src_det, .true.)
        write(6,*) 'Expecting ', nexcit, "excitations"
        call GenExcitations3 (src_det, ilut, det, flag, ex, par, found_all, &
                              .false.)
        do while (.not. found_all)
            ndet = ndet + 1
            call EncodeBitDet (det, det_list(:,ndet))

            call GenExcitations3 (src_det, ilut, det, flag, ex, par, &
                                  found_all, .false.)
        end do
        if (ndet /= nexcit) &
            call stop_all(this_routine,"Incorrect number of excitations found")

        ! Sort the dets, so they are easy to find by binary searching
        call sort(det_list, ilut_lt, ilut_gt)

        ! Lists to keep track of things
        allocate(generated_list(nexcit))
        allocate(contrib_list(nexcit))
        generated_list = .false.
        contrib_list = 0

        ! Repeated generation, and summing-in loop
        psingles = 0.0
        pdoubles = 1.0
        ngen = 0
        contrib = 0
        do i = 1, iterations
            if (mod(i, 10000) == 0) &
                write(6,*) i, '/', iterations, ' - ', contrib / real(ndet*i,dp)

            call gen_excit_4ind_weighted2 (src_det, ilut, det, tgt_ilut, 2, &
                                           ic, ex, par, pgen, helgen, store)
            if (det(1) == 0) cycle

            call EncodeBitDet (det, tgt_ilut)
            pos = binary_search(det_list, tgt_ilut, NIfD+1)
            if (pos < 0) then
                write(6,*) det
                write(6,'(b64)') tgt_ilut(0)
                write(6,*) 'FAILED DET', tgt_ilut
                call writebitdet(6, tgt_ilut, .true.)
                call stop_all(this_routine, 'Unexpected determinant generated')
            else
                generated_list(pos) = .true.

                ! Count this det, and sum in its contribution.
                ngen = ngen + 1
                contrib = contrib + 1.0_dp / pgen
                contrib_list(pos) = contrib_list(pos) + 1.0_dp / pgen
            end if
        end do

        ! How many of the iterations generated a good det?
        write(6,*) ngen, " dets generated in ", iterations, " iterations."
        write(6,*) 100_dp * (iterations - ngen) / real(iterations,dp), &
                   '% abortion rate'
        ! Contribution averages
        write(6, '("Averaged contribution: ", f15.10)') &
                contrib / real(ndet * iterations,dp)

        ! Output the determinant specific contributions
        iunit = get_free_unit()
        open(iunit, file="contribs_4ind", status='unknown')
        do i = 1, ndet
            call writebitdet(iunit, det_list(:,i), .false.)
            write(iunit, *) contrib_list(i) / real(iterations, dp)
        end do
        close(iunit)

        ! Check that all of the determinants were generated!!!
        if (.not. all(generated_list)) then
            write(6,*) count(.not.generated_list), '/', size(generated_list), &
                       'not generated'
            found_all = .true.
            do i = 1, ndet
                if (.not. generated_list(i)) then
                    call decode_bit_det(det, det_list(:,i))
                    hel = get_helement(src_det, det, ilut, det_list(:,i))
                    if (abs(hel) > 1.0e-6) then
                        found_all = .false.
                        call writebitdet(6, det_list(:,i), .false.)
                        write(6,*) hel
                    end if
                end if
            end do
            if (.not. found_all) &
                call stop_all(this_routine, "Determinant not generated")
        end if
        if (any(abs(contrib_list / iterations - 1.0_dp) > 0.01_dp)) then
            do i = 1, ndet
                call writebitdet(6, det_list(:,i), .false.)
                write(6,*) contrib_list(i) / (iterations - 1.0_dp)
            end do
            call stop_all(this_routine, "Insufficiently uniform generation")
        end if

        ! Clean up
        deallocate(det_list)
        deallocate(contrib_list)
        deallocate(generated_list)

    end subroutine


end module
