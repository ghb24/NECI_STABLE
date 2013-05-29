! Copyright (c) 2013, Ali Alavi
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
#include "macros.h"

! A new implementation file for csfs
module csf
    use constants, only: sizeof_int
    use systemdata, only: nel, brr, ecore, alat, nmsh, nbasismax, G1, nbasis,&
                          LMS, iSpinSkip, STOT, ECore, modk_offdiag
    use memorymanager, only: LogMemAlloc, LogMemDealloc
    use integralsdata, only: umat, fck, nmax
    use constants, only: dp, n_int, lenof_sign
    use dSFMT_interface, only: genrand_real2_dSFMT
    use sltcnd_mod, only: sltcnd, sltcnd_2
    use DetBitOps, only: EncodeBitDet, FindBitExcitLevel, count_open_orbs, &
                         get_bit_open_unique_ind, FindSpatialBitExcitLevel, &
                         DetBitEq
    use CalcData, only: InitiatorWalkNo
    use OneEInts, only: GetTMatEl
    use Integrals_neci, only: GetUMatEl
    use UMatCache, only: gtID
    use csf_data
    use timing_neci
    use util_mod, only: swap, choose
    use bit_rep_data, only: extract_sign
    use bit_reps, only: NIfD, NIfTot, NIfY

    implicit none

contains

    integer function num_csf_dets (ncsf)

        ! Return the number of determinants given a specified number of csfs.
        ! Note that if nopen == 0, ncsf == 0, but there is still a closed
        ! shell determinant that legitimately exists.
        !
        ! In: ncsf - The number of csfs (as from get_num_csfs)

        integer, intent(in) :: ncsf

        if (ncsf == 0) then
            num_csf_dets = 1
        else
            num_csf_dets = ncsf
        endif
    end function

    function CSFGetHelement (nI, nJ) result(hel_ret)
        
        ! Calculate the H-matrix element between two CSFs (nI, nJ)
        ! This is a wrapper function to allow working arrays to be on the
        ! stack.
        ! This is a wrapper for where we don't know ilutI,J; or for old code.
        !
        ! In:  nI, nJ   - The determinants to consider
        ! Ret: hel_ret  - The H-matrix element

        integer, intent(in) :: NI(nel), NJ(nel)
        integer(n_int) :: iLutI(0:NIfTot), iLutJ(0:NIfTot)
        HElement_t :: hel_ret

        ! These are needed just for the interface with get_csf_helement
        integer :: ic, ex(2,2)
        logical :: tParity
        HElement_t :: HElGen

        call EncodeBitDet (nI, iLutI)
        call EncodeBitDet (nJ, iLutJ)
        hel_ret = get_csf_helement(nI, nJ, iLutI, iLutJ, ic, ex, tParity, &
                                   HElGen)
    end function CSFGetHelement

    function get_csf_helement (nI, nJ, iLutI, iLutJ, notic, notex, &
                               nottParity, notHElGen) result (hel_ret)
        integer, intent(in) :: nI(nel), nJ(nel)
        integer, intent(in) :: notic, notex(2,2)
        integer(kind=n_int), intent(in) :: iLutI(0:NIfTot), iLutJ(0:NIfTot)
        logical, intent(in) :: nottParity
        HElement_t :: hel_ret
        HElement_t , intent(in) :: notHElGen

        integer :: nopen(2), nclosed(2), nup(2), ndets(2), IC, i
        integer :: S(2), Ms(2), iUnused
        logical :: bCSF(2), bBothCSF, tUnused

        ! Avoid compiler warnings
        iUnused=notic; iUnused=notex(1,1); tUnused=nottParity

        ! Are these both CSFs?
        bCSF(1) = iscsf(nI)
        bCSF(2) = iscsf(nJ)
        bBothCSF = bCSF(1) .and. bCSF(2)

        if ( (.not. bCSF(1)) .and. (.not. bCSF(2)) ) then
            ! Once again, pass things through better
            hel_ret = sltcnd (nI, iLutI, iLutJ, IC)
            if (IC == 0) then
                hel_ret = hel_ret + ECore
            else if (modk_offdiag) then
                hel_ret = -abs(hel_ret)
            end if
            return
        endif

        ! If the CSFs differ by more than 2 spin orbitals, the =0
        if (.not. bBothCSF) then
            IC = FindSpatialBitExcitLevel (iLutI, iLutJ)
        else
            IC = FindBitExcitLevel (iLutI, iLutJ)
        endif

        if (IC > 2) then
            hel_ret = (0)
            return
        endif

        ! Obtain statistics for each of the CSFs.
        S = 0
        if (bCSF(1)) &
            call get_csf_data(NI, ilutI, nopen(1), nclosed(1), S(1), Ms(1))
        if (bCSF(2)) &
            call get_csf_data(NJ, ilutJ, nopen(2), nclosed(2), S(2), Ms(2))

        ! If S are not consistent, then return 0
        ! Assume Ms is maintained as Ms==S
        ! .or. (Ms(1).ne.Ms(2))) then
        if (bBothCSF .and. (S(1).ne.S(2))) then
            hel_ret = (0) 
            return
        endif

        ! Use the maximal Ms value that we can (fewest determinants 
        ! required)
        Ms(1) = minval(S)
        Ms(2) = Ms(1)

        ! Get electronic details
        ! Using S instead of Ms to calculate nup, as this has the fewest
        ! determinants, and the Ms=S case is degenerate.
        nup = (nopen + S)/2
        do i=1,2
            if (bCSF(i)) then
                ndets(i) = int(choose(nopen(i), nup(i)))
            else
                ndets(i) = 1
            endif
        enddo

        ! Perform the calculation
        if (bBothCSF) then
            hel_ret = get_csf_helement_local (nI, nJ, iLutI, iLutJ, nopen, &
                                              nclosed, nup, ndets, IC)
        else if (bCSF(1)) then
            hel_ret = get_csf_helement_det (nI, iLutJ, nopen, nclosed, nup, &
                                            ndets)
        else
            call swap(nopen(1), nopen(2))
            call swap(nclosed(1), nclosed(2))
            call swap(nup(1), nup(2))
            call swap(ndets(1), ndets(2))
            hel_ret = get_csf_helement_det (nJ, iLutI, nopen, nclosed, nup, &
                                            ndets)
        endif

        if (modk_offdiag) then
            if (IC /= 0 .or. .not. DetBitEq(ilutI, ilutJ)) then
                hel_ret = -abs(hel_ret)
            end if
        end if

    end function

    function get_csf_helement_local (nI, nJ, iLutI, iLutJ, nopen, nclosed, &
                                     nup, ndets, IC) &
                                     result(hel_ret)

        ! The worker function for the above wrapper for calculating the
        ! H-matrix elements between two CSFs. By using a wrapper in this way,
        ! we can easily place the working arrays on the stack rather than
        ! using heap allocation (and therefore logging), or using pure
        ! function spaghetti code in the variable declarations.

        integer, intent(in) :: nI(nel), nJ(nel), nopen(2), nclosed(2)
        integer(kind=n_int), intent(in) :: iLutI(0:NIfTot), iLutJ(0:NIfTot)
        integer, intent(in) :: nup(2), ndets(2), IC
        HElement_t :: hel_ret

        ! Working arrays. Sizes calculated in calling function.
        integer :: yama1(nopen(1)), yama2(nopen(2))
        real(dp) :: coeffs1(ndets(1)), coeffs2(ndets(2))
        integer :: dets1(nel, ndets(1)), dets2(nel, ndets(2))
        integer(kind=n_int) :: ilut1(0:NIfTot,ndets(1)), ilut2(0:NIfTot,ndets(2))
        integer :: det_sum(ndets(1))

        integer :: det, i, j
        HElement_t :: sum1, Hel
        type(timer), save :: hel_timer0, hel_timer2, hel_timer1, hel_timer

        hel_timer%timer_name = 'Hel_timer'
        hel_timer0%timer_name = 'hel_timer0'
        hel_timer1%timer_name = 'hel_timer1'
        hel_timer2%timer_name = 'hel_timer2'

        call set_timer (hel_timer)

        ! Extract the Yamanouchi symbols from the CSFs
        call get_csf_yama (nI, yama1, nopen(1))
        call get_csf_yama (nJ, yama2, nopen(2))

        ! Depending on the number of spatial orbitals we differ by, call
        ! different optimised routines.
        if (IC == 2) then
            call halt_timer(hel_timer)
            call set_timer(hel_timer2)
            hel_ret = get_csf_helement_2 (nI, nJ, iLutI, iLutJ, nopen,  &
                                    nclosed, nup, ndets, dets1, dets2,&
                                    yama1, yama2, coeffs1, coeffs2)
            call halt_timer(hel_timer2)
            return
        endif

        ! The GENERAL routine follows now. Understand this before attempting
        ! to understand the optimised routines.

        ! Calculate all possible permutations to construct determinants
        ! (Where 0=alpha, 1=beta when generating NI, NJ below)
        if (IC == 0) then
            call csf_get_dets (nopen(1), nup(1), ndets(1), nel, dets1,det_sum)
        else
            call csf_get_dets (nopen(1), nup(1), ndets(1), nel, dets1)
        endif
        if ((nopen(1).eq.nopen(2)) .and. (nup(1).eq.nup(2))) then
            dets2 = dets1
        else
            call csf_get_dets (nopen(2), nup(2), ndets(2), nel, dets2)
        endif

        ! Get the coefficients
        do det=1,ndets(1)
            coeffs1(det) = csf_coeff(yama1,dets1(nclosed(1)+1:nel,det),&
                                     nopen(1))
        enddo
        do det=1,ndets(2)
            coeffs2(det) = csf_coeff(yama2,dets2(nclosed(2)+1:nel,det),&
                                     nopen(2))
        enddo

        ! If IC==0, then use the optimised routine for that case.
        call halt_timer(hel_timer)
        if (IC == 0) then
            call set_timer (hel_timer0)
            hel_ret = get_csf_helement_0 (nI, nopen(1), nclosed(1), nup(1), &
                                          ndets(1), coeffs1, coeffs2, dets1,&
                                          det_sum, all(yama1 == yama2))
            call halt_timer (hel_timer0)
            return
        endif

        call set_timer (hel_timer1)
        ! Generate determinants from spatial orbitals specified in NI, NJ
        ! We should be able to optimise this further.
        do det = 1,ndets(1)
            dets1(1:nclosed(1),det) = iand(NI(1:nclosed(1)), csf_orbital_mask)
            dets1(nclosed(1)+1:nel,det) = &
                    csf_alpha_beta(NI(nclosed(1)+1:nel), &
                                   dets1(nclosed(1)+1:nel,det))
            call EncodeBitDet (dets1(:,det), ilut1(:,det))
        enddo
        do det = 1,ndets(2)
            dets2(1:nclosed(2),det) = iand(NJ(1:nclosed(2)), &
                                           csf_orbital_mask)
            dets2(nclosed(2)+1:nel,det) = &
                    csf_alpha_beta(NJ(nclosed(2)+1:nel),&
                                   dets2(nclosed(2)+1:nel,det))
            call EncodeBitDet (dets2(:,det), ilut2(:,det))
        enddo

        ! Sum in all of the H-matrix terms between each of the component
        ! determinants.
        hel_ret = (0)
        do i=1,ndets(1)
            if (coeffs1(i) /= 0) then
                sum1 = (0)
                do j=1,ndets(2)
                    if (coeffs2(j) /= 0) then
                        Hel = sltcnd (dets1(:,i), ilut1(:,i), ilut2(:,j))
                        sum1 = sum1 + Hel * (coeffs2(j))
                    endif
                enddo
                hel_ret = hel_ret + sum1*(coeffs1(i))
            endif
        enddo

        call halt_timer(hel_timer1)
    end function

    function get_csf_helement_det (nI, iLutJ, nopen, nclosed, &
                                   nup, ndets) &
                                   result(hel_ret)

        ! The worker function for the above wrapper for calculating the
        ! H-matrix elements between a CSF and a normal determinant. By using 
        ! a wrapper in this way, we can easily place the working arrays on the
        ! stack rather than using heap allocation (and therefore logging), or 
        ! using pure function spaghetti code in the variable declarations.
        !
        ! n.b. nI/iLutI indicate the CSF, nJ/iLutJ --> The determinant.
        !      We assume that Ms/S etc. are correct.
        ! TODO: make nopen(2) --> nopen etc.

        integer, intent(in) :: nI(nel), nopen(2), nclosed(2)
        integer(kind=n_int), intent(in) :: iLutJ(0:NIfTot)
        integer, intent(in) :: nup(2), ndets(2)
        HElement_t :: hel_ret

        ! Working arrays. Sizes calculated in calling function.
        integer :: yama(nopen(1))
        real(dp) :: coeffs(ndets(1))
        integer :: dets(nel, ndets(1))
        integer(kind=n_int) :: ilut(0:NIfTot,ndets(1))

        integer :: det, i
        HElement_t :: hel

        ! Extract the Yamanouchi symbol from the CSF
        call get_csf_yama (nI, yama, nopen(1))

        ! Calculate all possibel permutations to construct determinants
        ! (Where 0=alpha, 1=beta below)
        call csf_get_dets (nopen(1), nup(1), ndets(1), nel, dets)

        ! Get the coefficients
        do det=1,ndets(1)
            coeffs(det) = csf_coeff (yama, dets(nclosed(1)+1:nel, det), &
                                     nopen(1))
        enddo

        ! Generate determinants from spatial orbitals specified in nI
        ! TODO: Can we optimise this further?
        do det=1,ndets(1)
            dets(1:nclosed(1),det) = iand(nI(1:nclosed(1)), csf_orbital_mask)
            dets(nclosed(1)+1:nel,det) = &
                   csf_alpha_beta(nI(nclosed(1)+1:nel), &
                                  dets(nclosed(1)+1:nel,det))
            call EncodeBitDet (dets(:,det), ilut(:,det))
        enddo

        ! Sum in all of the H-matrix terms
        hel_ret = (0)
        do i=1,ndets(1)
            if (coeffs(i) /= 0) then
                hel = sltcnd (dets(:,i), ilut(:,i), ilutJ)
                hel_ret = hel_ret + hel*(coeffs(i))
            endif
        enddo

    end function

    function get_csf_helement_0 (nI, nopen, nclosed, nup, ndets, coeffs1, &
                                 coeffs2, dets1, det_sum, bEqual) &
                                 result(hel_ret)
        
        ! The local worker function for calculating Helements between two CSFs
        ! which differ by 0 spatial orbitals (i.e. only differ by Yamanouchi
        ! symbols, or not at all).

        integer, intent(in) :: nI(nel)
        integer, intent(in) :: nopen, nclosed, nup, ndets
        real(dp), intent(in) :: coeffs1(ndets), coeffs2(ndets)
        integer, intent(in) :: dets1(nel,ndets)
        integer, intent(in) :: det_sum(ndets)
        logical, intent(in) :: bEqual
        HElement_t :: hel_ret

        integer :: nK(nel), id(nel), ex(2,2), elecs(2)
        integer :: ndown, idX, idN, ids, det, indj, i, j
        real(dp) :: diag_coeff
        HElement_t :: hel, hel2

        ! TODO: bEqual
        ! TODO: commenting

        ! Obtain spatial rather than spin indices
        nK = iand(nI, csf_orbital_mask)
        id = gtID(nK)

        ! Sum the coefficients of the diagonal matrix elements
        diag_coeff = sum(coeffs1*coeffs2)

        hel_ret = (0)
        if (diag_coeff /= 0) then
            ! Sum in the one electron integrals
            hel_ret = (diag_coeff) * &
                          sum(gettmatel(nK, nK))

            ! Sum in the terms which involve closed shell orbitals, which
            ! are the same across all involved determinants.
            hel = (0)
            do i=1,nclosed-1,2
                ! Within an orbital pair
                hel = hel + GetUmatEl(id(i), id(i), id(i), id(i))

                ! Between closed electron pairs
                if (i < nclosed-2) then
                    do j=i+2,nclosed-1,2
                        hel = hel + (4) * &
                                GetUmatEl(id(i), id(j), id(i), id(j))
                        hel = hel - (2) * &
                                GetUMatEl(id(i), id(j), id(j), id(i))
                    enddo
                endif

                ! Between a closed e- pair, and open electrons
                do j=nclosed+1,nel
                    idX = max(id(i), id(j))
                    idN = min(id(i), id(j))
                    hel = hel + (2) * &
                                getUMatEl(idN, idX,idN,idX)
                    hel = hel - GetUmatEl(idN, idX,idX,idN)
                enddo
            enddo
            hel_ret = hel_ret + (diag_coeff)*hel
        endif

        ! Sum in terms between orbitals pairs that differ
        hel = (0)
        do i=nclosed+1,nel-1
            do j=i+1,nel
                idX = max(id(i), id(j))
                idN = min(id(i), id(j))

                hel2 = GetUMatEl(idN, idX, idX, idN)
                
                do det=1,ndets
                    if (coeffs1(det) /= 0 .and. coeffs2(det) /= 0) then
                        ! Only include terms with matching Ms values.
                        ids = ieor(dets1(i,det), dets1(j,det))
                        if (ids == 0) then
                            hel = hel - hel2 * &
                                       (coeffs1(det)*coeffs2(det))
                        endif
                    endif
                enddo

                if (diag_coeff /= 0) then
                    hel = hel + (diag_coeff) * &
                                GetUMatEl(idN, idX, idN, idX)
                endif
            enddo
        enddo
        hel_ret = hel_ret + hel

        ! Now consider the cross terms (between differing component dets).
        ! Only cases differing by two spin orbitals count --> generate these.
        ! To maintain Ms and the spatial arrangement, can only swap two
        ! electron pairs with opposite Ms (i.e. AB --> BA)
        ndown = nel - nup
        do det=1,ndets
            elecs = -1
            call det_hel_2_pair (elecs(1), elecs(2), nopen, &
                                 dets1(nclosed+1:nel,det))

            do while (elecs(1) /= -1)
                ! The slow bit. Which determinant is this (to index coeffs
                ! array). 
                ! TODO: Is it quicker to re-calc the coeff?

                ! Get the position in list of dets. Alternative to det_pos.
                indj = det_perm_pos (dets1(nclosed+1:nel,det), elecs, &
                                     det_sum, nopen, ndets, det)

                if ( (coeffs1(det) /= 0 .and. coeffs2(indj) /= 0) .or. &
                     (coeffs1(indj) /= 0 .and. coeffs2(det) /= 0) ) then

                    ! Generate the excitation matrix.
                    do j=1,2
                        if (dets1(nclosed+elecs(j),det) == 1) then
                            ex(1,j) = get_beta(nK(nclosed + elecs(j)))
                        else
                            ex(1,j) = get_alpha(nK(nclosed + elecs(j)))
                        endif
                    enddo
                    ex(2,:) = ab_pair(ex(1,:))

                    ! We know this is a double spin-orbital excitation.
                    hel = sltcnd_2 (ex, .false.)

                    hel_ret = hel_ret &
                              + hel*(coeffs1(det)*coeffs2(indj))&
                              + hel*(coeffs1(indj)*coeffs2(det))
                endif

                call det_hel_2_pair (elecs(1), elecs(2), nopen, &
                                     dets1(nclosed+1:nel,det))
            enddo
        enddo

        ! If this a diagonal matrix element, sum in ECore
        if (bEqual) then
            hel_ret = hel_ret + (ECore)
        endif
        
    end function

    function get_csf_helement_2 (nI, nJ, iLutI, iLutJ, nopen, nclosed, &
                                 nup, ndets, dets1, dets2, yama1, yama2, &
                                 coeffs1, coeffs2) result(hel_ret)

        ! Given a case where IC == 2, calculate the . Thus we can make
        ! some rather stark approximations. In particular, as we definitely
        ! differ by 2 spatial oribitals, only those orbitals can differ in Ms
        ! value, otherwise the overall pair differs my more than 2 spin
        ! orbitals --> Problematic.
        !
        ! We assign the fastest changing bits in the permutations (in
        ! csf_get_dets_ind) to the differing spatial orbitals. Thus all of the
        ! common parts occur in the same order for both CSFs (although if
        ! nopen varies, we may need to skip some extra terms in one or the
        ! other). This makes the matrix nearly diagonal, so we can sum over it
        ! as such --> only calculate terms which could possibly be non-zero.
        !
        ! TODO: shift declaratino of dets into here --> we only need it to be
        !       nopen long, as we don't generate the dets.

        use constants, only: bits_n_int
        integer, intent(in) :: nI(nel), nJ(nel), nopen(2), nclosed(2)
        integer(kind=n_int), intent(in) :: iLutI(0:NIfTot), iLutJ(0:NIfTot)
        integer, intent(in) :: nup(2), ndets(2)
        integer, intent(in) :: yama1(nopen(1)), yama2(nopen(2))
        real(dp), intent(out) :: coeffs1(ndets(1)), coeffs2(ndets(2))
        integer, intent(inout) :: dets1(nel, ndets(1)), dets2(nel, ndets(2))
        HElement_t :: hel_ret, hel, sum1, umatel(2)

        integer :: nop_uniq(2), det, i, j, k, l, m, ex(2,2), &
                   uniq_id(4,2), ms1(ndets(1)), ms2(ndets(2)), &
                   nsign(2), tsign_id(4,2), id(2,2), &
                   ex_ms_ind(2,2), ex_ms(2,2)

        logical :: dets_change1(ndets(1)), dets_change2(ndets(2)), tSign, &
                   delta_tsign1(ndets(1)), delta_tsign2(ndets(2)), tSign_tmp

        ! count the number of singles which differ between the two dets, and
        ! get their indices in the open section. Also, get the indices of any
        ! orbitals whose Ms value may affect tSign.
        call get_bit_open_unique_ind (iLutI, iLutJ, uniq_id, nop_uniq, &
                                      tsign_id, nsign, 2)

        ! Get the excitation matrix, and more importantly the parity of the
        ! 'standard' excitation (i.e. all singles are beta).
        ex(1,1) = 2
        call GetBitExcitation (iLutI, iLutJ, ex, tSign)

        ! Calculate all possible permutations to construct determinants
        ! (Where 0=alpha, 1=beta when generating NI, NJ below)
        call csf_get_dets_ind (nopen(1), nup(1), ndets(1), nel, &
                               nop_uniq(1), uniq_id(:,1), tsign_id(:,1), &
                               nsign(1), dets1, ms1, delta_tsign1)
        call csf_get_dets_ind (nopen(2), nup(2), ndets(2), nel, &
                               nop_uniq(2), uniq_id(:,2), tsign_id(:,2), &
                               nsign(2), dets2, ms2, delta_tsign2)

        ! TODO: can we do this in the previous bit?
        ! Mark each permutation for which each of the differing orbitals are
        ! permuted in a canonical order --> the end of the set which mixes.
        call mark_change_2 (dets_change1, dets1, nclosed(1), ndets(1), &
                            nop_uniq(1), uniq_id(:,1))
        call mark_change_2 (dets_change2, dets2, nclosed(2), ndets(2), &
                            nop_uniq(2), uniq_id(:,2))

        ! For each of the orbitals involved, are we exciting to/from a double,
        ! or a single. If appropriate, what is the index into the determinant
        ! too be able to examine the details.
        do i=1,2
            if (is_in_pair(ex(i,1), ex(i,2))) then
                ex_ms_ind(i,:) = 0
                ex_ms (i,1) = -1
                ex_ms (i,2) = 1
            else
                ! TODO: tidy up i/j separation, i.e. coalesce the loops.
                do j=1,2
                    if (IsOcc(iLutI, ab_pair(ex(i,j)))) then
                        do k=1,nop_uniq(3-i)
                            if (i == 1) then
                                if (is_in_pair(ex(i,j), iand(nJ(nclosed(2)+uniq_id(k,2)), csf_orbital_mask))) then
                                    ex_ms_ind(i,j) = -nclosed(2)-uniq_id(k,2)
                                    exit
                                endif
                            else
                                if (is_in_pair(ex(i,j), iand(nI(nclosed(1)+uniq_id(k,1)), csf_orbital_mask))) then
                                    ex_ms_ind(i,j) = nclosed(1)+uniq_id(k,1)
                                    exit
                                endif
                            endif
                        enddo
                    else
                        do k=1,nop_uniq(i)
                            if (i == 1) then
                                if (is_in_pair(ex(i,j), iand (nI(nclosed(1)+uniq_id(k,1)), csf_orbital_mask))) then
                                    ex_ms_ind(i,j) = nclosed(1) + uniq_id(k,1)
                                    exit
                                endif
                            else
                                if (is_in_pair(ex(i,j), iand (nJ(nclosed(2)+uniq_id(k,2)), csf_orbital_mask))) then
                                    ex_ms_ind(i,j) = -nclosed(2)-uniq_id(k,2)
                                    exit
                                endif
                            endif
                        enddo
                    endif
                enddo
            endif
        enddo


        ! Get the coefficients
        do det=1,ndets(1)
            coeffs1(det) = csf_coeff(yama1,dets1(nclosed(1)+1:nel,det),&
                                     nopen(1))
        enddo
        do det=1,ndets(2)
            coeffs2(det) = csf_coeff(yama2,dets2(nclosed(2)+1:nel,det),&
                                     nopen(2))
        enddo

        ! Obtain the two UMAT components which may get used in the sum
        id = gtID(ex)
        umatel(1) = GetUMATEl (id(1,1), id(1,2), id(2,1), id(2,2))
        umatel(2) = GetUMATEL (id(1,1), id(1,2), id(2,2), id(2,1))

        ! Loop through all of the terms where all but the fastest changing
        ! bits of each component determinant are the same. Slightly nasty
        ! looping code... --> the matrix is block diagonalised.
        j = 1
        i = 1
        hel_ret = (0)
        do while (i <= ndets(1) .and. j <= ndets(2))

            ! Avoid the extra terms which appear if the nopen values are
            ! different (more flexibility to maintain Ms in the permutations)
            if (ms1(i) /= ms2(j)) then
                if (nopen(1) > nopen(2)) then
                    do while (i <= ndets(1))
                        i = i + 1
                        if (dets_change1(i-1)) exit
                    enddo
                else
                    do while (j <= ndets(2))
                        j = j + 1
                        if (dets_change2(j-1)) exit
                    enddo
                endif
                cycle
            endif
            
            ! For each valid determinant in CSF1, loop through all of the
            ! connected determinants on the other determinant.
            if (coeffs1(i) /= 0) then
                sum1 = (0)
                k = j
                do while (k <= ndets(2))
                    if (coeffs2(k) /= 0) then

                        ! Obtain the Ms matrix for the excitation (i.e. the Ms
                        ! values of the orbitals excited from/to)
                        do l=1,2
                        do m=1,2
                            if (ex_ms_ind(l,m) < 0) then
                                ex_ms(l,m) = 2*dets2(abs(ex_ms_ind(l,m)),k)-1
                                if (l == 1) ex_ms(l,m) = -ex_ms(l,m)
                            else if (ex_ms_ind(l,m) > 0) then
                                ex_ms(l,m) = 2*dets1(ex_ms_ind(l,m),i)-1
                                if (l == 2) ex_ms(l,m) = -ex_ms(l,m)
                            endif
                        enddo
                        enddo

                        ! Include only those terms allowed by Ms symmetry.
                        hel = (0)
                        if ( (ex_ms(1,1) == ex_ms(2,1)) .and. &
                             (ex_ms(1,2) == ex_ms(2,2)) ) then
                            hel = hel + umatel(1)
                        endif
                        if ( (ex_ms(1,1) == ex_ms(2,2)) .and. &
                             (ex_ms(1,2) == ex_ms(2,1)) ) then
                            hel = hel - umatel(2)
                        endif

                        ! Calculate the parity of this pair of dets.
                        tSign_tmp = tSign .neqv. (delta_tsign1(i) .neqv. &
                                                 delta_tsign2(k))
                        if (tSign_tmp) hel = -hel

                        ! Update the overall sum.
                        sum1 = sum1 + hel*(coeffs2(k))
                    endif

                    ! Increment the second index, but remain within the block
                    k = k + 1
                    if (dets_change2(k-1)) exit
                enddo
                hel_ret = hel_ret + sum1*(coeffs1(i))
            endif

            ! If appropriate, move the second index on to the next block
            ! before incrementing the first index (i.e. progress both blocks
            ! to the next one at the same time).
            if (dets_change1(i)) then
                do while (j <= ndets(2))
                    j = j + 1
                    if (dets_change2(j-1)) exit
                enddo
            endif

            i = i + 1
        enddo

    end function

    subroutine det_hel_2_pair (elecA, elecB, nopen, det)
        
        ! Pick pairs of electrons from the specified 'determinant'. The
        ! determinant contains only 1s and 0s, relating to the alpha and beta
        ! electrons. Pick pairs only in a standard order (1s before 0s),
        ! avoids any double counting in get_helement_0
        !
        ! In:    nopen        - Number of open electrons to consider
        !        dets         - Array of 1s and 0s
        ! InOut: elecA, elecB - Contains the electron pair to consider.
        !                       -1 initialises, and indicates the last pair

        integer, intent(in) :: nopen, det(nopen)
        integer, intent(inout) :: elecA, elecB

        ! Are we initialising
        if (elecA == -1) then
            do elecA=1,nopen
                if (det(elecA) == 1) exit
            enddo

            ! There are no electrons in the A position.
            if (elecA > nopen) then
                elecA = -1
                elecB = -1
            endif
        else
            ! Increment the B electron
            do elecB=elecB+1,nopen
                if (det(elecB) == 0) exit
            enddo

            ! If we have run out of B electrons, increment the A electron, and
            ! signal that we need to start again with the B elecs.
            if (elecB > nopen) then
                elecB = -1
                do elecA=elecA+1,nopen
                    if (det(elecA) == 1) exit
                enddo

                ! Have we used up all possibilities?
                if (elecA > nopen) elecA = -1
            endif
        endif

        ! Do we need to find another B electron?
        if (elecB == -1 .and. elecA /= -1) then
            do elecB=elecA,nopen
                if (det(elecB) == 0) exit
            enddo

            ! There are no possbile B electrons for this A (or any later ones)
            ! --> end of string.
            if (elecB > nopen) then
                elecA = -1
                elecB = -1
            endif
        endif
    end subroutine

    subroutine mark_change_2 (dets_change, dets, nclosed, ndets, nop_uniq,&
                              uniq_id)

        ! Look through the specified determinants, after they have been filled
        ! with all possible permutations. If the most rapidly varying bits
        ! (nop_uniq of them, specified in uniq_id) are in a canonical order
        ! (all 1s then 0s), then mark that position in dets_change as true.
        ! Helper function for get_csf_helement_2.
        ! 
        ! In:  dets        - The array of determinants, open orbitals filled
        !                    with permutations summing to 2Ms.
        !      nclosed     - Number of paired electrons
        !      ndets       - Number of determinants
        !      nop_uniq  - Number of unique singles --> no of rapidly
        !                    varying terms to consider
        !      uniq_id       - Indices of rapidly varying terms - nclosed.
        ! Out: dets_change - true if det in canonical order, otherwise false.

        integer, intent(in) :: nclosed, ndets, nop_uniq
        integer, intent(in) :: dets(nel, ndets), uniq_id(4)
        logical, intent(out) :: dets_change(ndets)
        logical :: bChange
        integer :: i, j

        dets_change = .true.
        if (nop_uniq == 0) return

        do i=1,ndets
            bChange = .false.
            do j=1,nop_uniq
                if (.not. bChange .and. dets(uniq_id(j)+nclosed,i) == 0) &
                    bChange = .true.

                if (bChange .and. dets(uniq_id(j)+nclosed,i) /= 0) then
                    dets_change(i) = .false.
                    exit
                endif
            enddo                    
        enddo

    end subroutine

    subroutine csf_get_bit_perm (nopen, nup, ndets, ilut)
    
        ! As for csf_get_dets, but working with a bit representation of the
        ! permutation rather than separate integers
        
        use constants, only: bits_n_int
        integer, intent(in) :: ndets, nup, nopen
        integer(kind=n_int), intent(out) :: ilut(NIfY,ndets)
        integer :: i, det, comb(nup)
        logical :: bInc

        if (nopen == 0) return

        forall (i=1:nup) comb(i) = i
        ilut = 1
        do det = 1, ndets
            forall (i=1:nup) ilut((comb(i)-1)/bits_n_int,det) = &
                            ibclr(ilut((comb(i)-1)/bits_n_int,det), mod(comb(i)-1,bits_n_int))

            do i=1,nup
                bInc = .false.
                if (i == nup) then
                    bInc = .true.
                else if (i < nup) then
                    if (comb(i+1) /= comb(i) + 1) bInc = .true.
                endif

                if (bInc) then
                    comb(i) = comb(i) + 1
                    exit
                else
                    comb(i) = i
                endif
            enddo
        enddo
    end subroutine

    integer function det_pos (det, nopen, ndown, perm)
        
        ! Obtain an index for the number of the permutation generated
        ! according to the ordering used in csf_get_dets.
        ! See: The art of Computer Programming, volume 4, Fascicle 3, pg. 
        !
        ! In:  det   - An array of 0s, 1s which have been permuted
        !      nopen - Size of array det
        !      ndown - Number of 0s
        !      perm  - Treat as if indices perm(1), perm(2) are swapped.
        ! Ret:       - Std. index of determinant.

        integer, intent(in) :: nopen, ndown
        integer, intent(in) :: det(nopen)
        integer, intent(in), optional :: perm(2)
        integer :: i, pos

        ! Start at 1 not 0, due to fortrans indexing...
        det_pos = 1
        pos = 0
        do i=1,ndown
            do pos=pos+1,nopen
                if (present(perm)) then
                    if (pos == perm(1)) then
                        if (det(perm(2)) == 0) exit
                    else if (pos == perm(2)) then
                        if (det(perm(1)) == 0) exit
                    else if (det(pos) == 0) then
                        exit
                    endif
                else if (det(pos) == 0) then
                    exit
                endif
            enddo
            if (pos > nopen) exit
        !    write (6, '(i5)', advance='no') pos

            det_pos = det_pos + int(choose(pos-1, i),sizeof_int)
        enddo
    end function

    function det_perm_pos (det, perm, det_sum, nopen, ndets, deti) result(pos)

        ! Given a permutation det, at position deti in the array of
        ! determinants generated by csf_get_dets, find the index of the
        ! permutation generated if we switch the values at indices perm(1:2)
        ! Note that det(perm(1)) /= det(perm(2)) (one must equal 1, the other
        ! must equal 2).
        ! 
        ! Currently this assumes that nopen <= bits_n_int. This is reasonable as a CSF
        ! of that size would be unfeasible to use - we assume that on such a
        ! large system we would be using TRUNCATE-CSF. This might change with
        ! use of the initiator algorithm, in which case we will have to
        ! generalise the code to use an array (but it would be slightly
        ! slower).
        !
        ! nb. This routine does not need to know nup/ndown
        ! nb. This acts as a faster alternative to using det_pos, which is a
        !     more general function.
        !
        ! In:  det     - The permutation of 1s/0s
        !      perm    - The indices to swap
        !      det_sum - The array of sum(2**(i-1)) forall i where index of 0
        !                in the determinant is 0. This monotonically increases
        !                through the set of determinants --> allows binary
        !                searching.
        !      nopen   - Number of unpaired electrons (size of det)
        !      ndets   - The possible number of permutations (size of det_sum)
        !      deti    - Index of det in det_sum
        ! Ret: pos     - New index in det_sum

    ! TODO: comment

        use constants, only: bits_n_int
        integer, intent(in) :: ndets, nopen, deti
        integer, intent(in) :: det(nopen), perm(2)
        integer, intent(in) :: det_sum(ndets)
        integer :: pos

        integer :: sumdet, hi, lo

        if (nopen > bits_n_int) call stop_all ("det_perm_pos", "nopen too large. Need&
                                   & to move to multi-integer representation")

        ! Calculate the new sum value
        sumdet = det_sum(deti)
        if (det(perm(1)) == 1) then
            sumdet = ibset(ibclr(sumdet, perm(2)-1), perm(1)-1)
        else
            sumdet = ibset(ibclr(sumdet, perm(1)-1), perm(2)-1)
        endif

        ! Prepare binary search range
        if (sumdet > det_sum(deti)) then
            hi = ndets
            lo = deti + 1
        else
            hi = deti - 1
            lo = 1
        endif

        ! Perform a binary search to find the item.
        do while (hi /= lo)
            pos = int(real(hi + lo) / 2)

            if (sumdet == det_sum(pos)) then
                exit
            else if (sumdet > det_sum(pos)) then
                lo = pos + 1
            else
                hi = pos - 1
            endif
        enddo

        if (hi == lo) pos = hi
    end function


    subroutine csf_get_dets (nopen, nup, ndets, nel, dets, pos_sum)

        ! Fill the last nopen electrons of each determinant with 0 (alpha) or
        ! 1 (beta) in all possible permutations with nup alpha electrons.
        !
        ! Out:  pos_sum - Sum of 2^(pos-1) for the locations of the zeros.
        !                 This monotonically increases through the set, and
        !                 thus allows binary searching.


        integer, intent(in) :: ndets, nup, nopen, nel
        integer, intent(out) :: dets (nel, ndets)
        integer, intent(out), optional :: pos_sum (ndets)
        integer comb(nup), i, j
        logical bInc

        if (nopen == 0) return

        forall (i=1:nup) comb(i) = i
        dets(nel-nopen+1:,:) = 1
        do i=1,ndets
            forall (j=1:nup) dets(nel-nopen+comb(j), i) = 0
            if (present(pos_sum)) then
                pos_sum(i) = sum(ibset(0, comb-1))
            endif
            do j=1,nup
                bInc = .false.
                if (j == nup) then
                    bInc = .true.
                else if (j < nup) then
                    if (comb(j+1) /= comb(j) + 1) bInc = .true.
                endif

                if (bInc) then
                    comb(j) = comb(j) + 1
                    exit
                else
                    comb(j) = j
                endif
            enddo
        enddo
    end subroutine

    subroutine csf_get_dets_ind (nopen, nup, ndets, nel, nop_uniq, uniq_id, &
                                 tsign_id, nsign, dets, ms, tSign)

        ! Fill the last nopen electrons of each determinant with 0 (alpha) or
        ! 1 (beta) in all possible permutations with nup alpha electrons.
        ! Specify a number of positions, which should be the fastest varying,
        ! in the array uniq_id.
        !
        ! Also, excluding the specified fast varying bits, count the number of
        ! bits which are set in the remainder of the permutation (equivalent
        ! to the sum of the Ms values of this region, even if not numerically
        ! equivalent).
        !
        ! In:  nopen    - Num. unpaired orbitals (size of permutations)
        !      nup      - Num. orbitals in 'alpha' state
        !      ndets    - Num. permutations to generate
        !      nel      - Num. electrons in total
        !      nop_uniq - Num. of orbitals to place as 'fastest'
        !      uniq_id  - Indices of fastest changing orbitals - nclosed
        !      tsign_id - Indices of terms which can affect parity of det.
        ! Out: dets  - Array, last nopen bits of each det contain permutations
        !      ms    - Sum of terms other than the specified fast varying ones
        !      tsign - Does the permutation change the parity (tSign) of the
        !              det from the 'standard' csf one (all betas)

        integer, intent(in) :: ndets, nup, nopen, nel, nop_uniq, uniq_id(*)
        integer, intent(in) :: tsign_id(*), nsign
        integer, intent(out) :: dets (nel,ndets), ms(ndets)
        logical, intent(out) :: tSign (ndets)
        integer comb(nup), i, j, id(nopen), pos, posu, nclosed
        logical bInc

        if (nopen.eq.0) then
            ms = 0
            tSign = .false.
            return
        endif

        nclosed = nel - nopen

        ! Generate mapping of positions to cause the unique indices to vary as
        ! fast as possible
        id(1:nop_uniq) = uniq_id(1:nop_uniq)
        pos = nop_uniq+1
        posu = 1
        do i=1,nopen
            if (posu > nop_uniq) exit
            if (i == uniq_id(posu)) then
                posu = posu + 1
            else
                id(pos) = i
                pos = pos + 1
            endif
        enddo
        if (i <= nopen) forall (j=i:nopen) id(j) = j

        ! Calculate permutations and place in dets
        forall (i=1:nup) comb(i) = i
        dets(nclosed+1:,:) = 1
        do i=1,ndets

            ! Write zeros to the correct place in the determinant.
            forall (j=1:nup) dets(nclosed + id(comb(j)), i) = 0

            ! Sum the terms in the fastest changing (nop_uniq) positions.
            ms(i) = nup - sum(dets(nclosed + id(1:nop_uniq), i))

            ! Should tSign vary from the 'standard' one? If an excitation
            ! from a double occurs from a beta orbital, it changes tSign (a
            ! permutation is required to get it into the std. CSF order). Thus
            ! the orbital left behind is alpha to get a change --> det==0.
            tsign(i) = btest(sum(dets(nclosed+tsign_id(1:nsign), i)),0)

            ! Generate the next permutation.
            do j=1,nup
                bInc = .false.
                if (j == nup) then
                    bInc = .true.
                else if (j < nup) then
                    if (comb(j+1) /= comb(j) + 1) bInc = .true.
                endif

                if (bInc) then
                    comb(j) = comb(j) + 1
                    exit
                else
                    comb(j) = j
                endif
            enddo
        enddo
    end subroutine

    pure subroutine csf_get_yamas (nopen, sfinal, yama, ncsf_max)
        
        ! Obtain all possible Yamanouchi symbols for for the system
        !
        ! In:  nopen    - Number of open shell electrons
        !      sfinal   - Desired final spin
        !      ncsf_max - Max number of csfs to generate (size of yama)
        ! Out: yama     - Array of Yamanouchi symbols

        integer, intent(in) :: sfinal
        integer, intent(in) :: nopen, ncsf_max
        integer, intent(out) :: yama (ncsf_max, nopen)

        integer spin (ncsf_max, nopen)
        integer npos, csf, ncsf, ncsf_next

        ! Empty Yamanouchi symbol of nopen == 0.
        if (nopen == 0 .or. ncsf_max == 0) return

        ! Walk through tree. Start at elec == nopen, and walk back through
        ! tree taking all possible routes (branch at every point where there
        ! is more than one valid route to the start.
        spin(1,nopen) = sfinal
        ncsf = 1
        ncsf_next = ncsf
        do npos = nopen, 2, -1
            do csf=1,ncsf
                if (spin(csf,npos) .lt. npos) then
                    spin(csf,npos-1) = spin(csf,npos) + 1
                    yama(csf,npos) = 2
                    if (spin(csf,npos) .ne. 0) then
                        ncsf_next = ncsf_next + 1
                        if (ncsf_next .gt. ncsf_max) exit
                        !spin(ncsf_next,npos:nopen) = spin(csf,npos:nopen)
                        yama(ncsf_next,npos+1:nopen) = yama(csf,npos+1:nopen)
                        spin(ncsf_next,npos-1) = spin(csf,npos) - 1
                        yama(ncsf_next,npos) = 1
                    endif
                else
                    spin(csf,npos-1) = spin(csf,npos) - 1
                    yama(csf,npos) = 1
                endif
            enddo
            ncsf = ncsf_next
            if (ncsf .gt. ncsf_max) exit
        enddo
        yama(:,1) = 1
    end subroutine
        
    integer elemental function csf_alpha_beta (num, det)

        ! Convert num (a member of CI) to an alpha or beta spin orbital where
        ! det==1 --> alpha

        integer, intent(in) :: num, det

        csf_alpha_beta = iand(num, csf_orbital_mask) - 1
        if (det .eq. 1) then
            csf_alpha_beta = ior(csf_alpha_beta, 1)
        else
!            csf_alpha_beta = iand(csf_alpha_beta,Z'fffffffe')
            csf_alpha_beta = iand(csf_alpha_beta,-2)
        endif
        csf_alpha_beta = csf_alpha_beta + 1
    end function

    subroutine get_csf_yama(nI, yama, nopen)

        ! Extract the Yamanouchi symbol from the supplied CSF
        ! This assumes that the passed yama array is the correct size, if
        ! it is too small, then the symbol will be truncated.

        integer, intent(in) :: nI(nel)
        integer, intent(in) :: nopen
        integer, intent(out) :: yama(nopen)
        integer i, nclosed

        nclosed = nel - nopen
        do i=1,nopen
            if (btest(NI(i + nclosed), csf_yama_bit)) then
                yama(i) = 1
            else
                yama(i) = 2
            endif
        enddo
    end subroutine

    subroutine get_csf_bit_yama (nI, yama)
        
        ! Extract the Yamanouchi symbol from the supplied CSF as part of
        ! a bit determinant - the same as part of the bit det obtained by
        ! EncodeBitDet.

        use constants, only: bits_n_int
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(out) :: yama(NIfY)
        integer :: i, pos, bit, nopen

        nopen = 0
        yama = 0
        do i=1,nel
            ! Only consider open electrons. The first has csf_yama_bit set.
            if (nopen > 0 .or. btest(nI(i), csf_yama_bit)) then
                pos = 1 + nopen / bits_n_int
                bit = mod(nopen, bits_n_int)

                if (btest(nI(i), csf_yama_bit)) &
                    yama(pos) = ibset(yama(pos), bit)
                nopen = nopen + 1
            endif
        enddo
    end subroutine

    subroutine csf_to_old_csf (nI, nJ)
        
        ! Convert the CSF to the old representation (Alex's rep.) for testing
        ! purposes.
        
        use legacy_data, only: CSF_NBSTART
        integer, intent(in) :: nI(nel)
        integer, intent(out) :: nJ(nel)
        integer :: i, nopen, nup, orb
        logical :: open_shell

        if (.not. iscsf(nI)) then
            nJ = nI
            return
        endif

        nup = 0
        nopen = 0
        do i=1,nel
            orb = iand(nI(i), csf_orbital_mask)

            ! Detect the start of the open shell region
            if (.not. open_shell .and. btest(nI(i),csf_yama_bit)) then
                open_shell = .true.
                nopen = nel - i + 1
                ! Aim for an Ms value of LMS
                nup = (nopen + LMS) / 2
                nopen = 0
            endif

            if (open_shell) then
                nopen = nopen + 1
                nJ(i) = 4*(orb - 1)
                if (btest(nI(i), csf_yama_bit)) nJ(i) = ibset(nJ(i),1)
                if (nopen <= nup) nJ(i) = ibset(nJ(i), 0)
                nJ(i) = nJ(i) + csf_nbstart
            else
                nJ(i) = orb
            endif
        enddo
    end subroutine

    subroutine get_csf_data(nI, ilut, nopen, nclosed, S, Ms)

        ! Obtains the number of open shell electrons, the total spin
        ! and the Ms value for the specified csf.
        ! NB. The Ms value is no longer being used. We ASSERT that Ms == STOT

        integer, intent(in) :: NI(nel) 
        integer(kind=n_int), intent(in) :: ilut(0:NIfTot)
        integer, intent(out) :: nopen, nclosed
        integer, intent(out) :: S, Ms
        integer i
        logical open_shell

        if (.not. bTest(nI(1), csf_test_bit)) then
            S = STOT
            Ms = LMS
            nopen = count_open_orbs(ilut)
            nclosed = nel - nopen
            return
        endif

        nopen = 0
        nclosed = 0
        S = 0
        Ms = 0
        open_shell = .false.
        do i=1,nel
            ! Closed shell until they are no longer paired.
            if ((.not.open_shell) .and. &
                btest(NI(i),csf_yama_bit)) open_shell = .true.
                
            if (.not. open_shell) then
                nclosed = nclosed + 1
            else
                if (btest(NI(i), csf_yama_bit)) then
                    S = S + 1
                else
                    S = S - 1
                endif
                if (btest(NI(i), csf_ms_bit)) then
                    Ms = Ms + 1
                else
                    Ms = Ms - 1
                endif
            endif
        enddo
        nopen = nel - nclosed
    end subroutine

    pure function get_num_csfs (nOpen, S) result (ncsf)

        ! Calculates the total number of CSFs possible for a system with
        ! nOpen unpaired electrons, and a total spin of S/2.
        ! This is the same as the number of available Serber functions
        !
        ! In:  nopen - Number of unpaired electrons
        !      S     - 2*total spin of system.
        ! Ret: ncsf  - Number of CSFs given spatial structure.

        integer, intent(in) :: nOpen
        integer, intent(in) :: S
        integer :: ncsf

        if ((nopen < 0) .or. (mod(nOpen+S, 2) /= 0))then
            ncsf = 0
        else
            ncsf = (2*S + 2) * int(choose(nOpen, (nOpen+S)/2),sizeof_int)
            ncsf = ncsf / (nOpen + S + 2)
        endif
    end function

    ! TODO: We can make this more efficient. Each random number should
    !       contain enough information to make more than 1 choice (split it
    !       up into a variety of bits).
    function random_spin_permute (spins, Ms) result (no_dets)

        ! Take the selection of spins (in whatever order they are given) and
        ! apply a random Ms value to them.
        !
        ! In:    Ms      - 2 * the required Ms value
        ! InOut: spins   - Array of (unpaired) spin orbitals
        ! Ret:   no_dets - The number of dets chosen from

        integer, intent(inout) :: spins(:)
        integer, intent(in) :: Ms
        integer :: no_dets 

        integer :: nopen, nup, nchoose, pos, i
        integer :: choice(ubound(spins,1)), perm(ubound(spins,1))
        real(dp) :: r

        ! How many alpha elecs do we have. If fewer than half, permute betas.
        nopen = ubound(spins, 1)
        nup = (nopen - Ms) / 2

        ! We want to make the fewest (random) selections possible.
        nchoose = min(nup, nopen - nup)

        forall (i=1:nopen) choice(i) = i

        ! Select nchoose positions at random, and place them at the end.
        do i = 1, nchoose
            r = genrand_real2_dSFMT()
            pos = int(real(nopen-i+1,dp)*r) + 1

            if (pos /= nopen-i+1) then
                call swap(choice(nopen-i+1), choice(pos))
            endif
        enddo

        ! Generate the permutation
        if (nchoose == nup) then
            perm(choice(1:nopen-nchoose)) = 0
            perm(choice(nopen-nchoose+1:nopen)) = 1
        else
            perm(choice(1:nopen-nchoose)) = 1
            perm(choice(nopen-nchoose+1:nopen)) = 0
        endif

        ! Generate the determinant, correctly sorted, with the specified
        ! alpha/beta structure.
        spins = iand(spins, csf_orbital_mask)
        spins = csf_alpha_beta(spins, perm)

        ! How many dets were there to choose from?
        no_dets = int(choose(nopen, nchoose),sizeof_int)

    end function

    ! TODO: We can make this more efficient. Each random number should
    !       contain enough information to make more than 1 choice (split it
    !       up into a variety of bits).
    integer function csf_get_random_det (nI, nopen, Ms)

        ! Convert the CSF or CSF ordered determinant nI into a normal
        ! determinant with spin orbitals which are any of the determinants
        ! which would be a component of the CSF.
        !
        ! In:    nopen - The number of open shell electrons
        !        Ms    - 2 * The required Ms value for the determinant
        ! InOut: nI    - The CSF to consider, and the determinant to return
        ! Ret:         - The number of determinants chosen from
        
        integer, intent(inout) :: nI(nel)
        integer, intent(in) :: nopen
        integer, intent(in) :: Ms

        csf_get_random_det =  random_spin_permute (nI(nel-nopen+1:nel), Ms)
        nI(1:nel-nopen) = iand(nI(1:nel-nopen), csf_orbital_mask)
        call csf_sort_det_block (nI, 1, nopen)

    end function

    function det_to_random_csf (nI) result(ncsf)

        ! From a determinant, generate a (random) CSF which has the same
        ! spatial structure
        !
        ! InOut: nI   - The determinant and output CSF
        ! Ret:   ncsf - The number of possible Yamanouchi symbols for this
        !               spatial structure

        integer, intent(inout) :: nI(nel)
        integer :: ncsf

        integer :: sings(nel)
        integer :: i, pos, nopen

        pos = 1
        nopen = 0
        i = 1
        do while (i <= nel .and. pos <= nel)

            ! Is this a double?
            if (i < nel .and. is_beta(nI(i))) then
                if (nI(i+1) == nI(i) + 1) then
                    nI(pos:pos+1) = nI(i:i+1)
                    i = i + 2
                    pos = pos + 2
                    cycle
                endif
            endif

            ! Otherwise it must be a single
            nopen = nopen + 1
            sings(nopen) = nI(i)
            i = i + 1
        enddo

        if (nopen == 0) then
            ncsf = 0
            return
        endif

        ! All unpaired electrons are stored as 'beta' in CSF representation.
        nI(nel-nopen+1:nel) = get_beta(sings(1:nopen))

        call csf_apply_random_yama (nI, nopen, STOT, ncsf, .false.)
    end function

    subroutine csf_apply_random_yama (nI, nopen, S, ncsf, tForceChange)

        ! Apply a random Yamanouchi symbol to the specified CSF or determinant
        ! with CSF ordering. Currently generates all possible Yamanouchi
        ! symbols and then picks.
        ! TODO: optimise, or store generated symbols (for trunc_csf)
        !
        ! In:  nI           - Integer representation of determinant/CSF.
        !      nopen        - Number of unpaired electrons.
        !      S            - The desired 2*total spin required.
        !      tForceChange - If ncsf > 1, ensure that we change the CSF.
        ! Out: ncsf         - Number of CSFs we have picked from.

        integer, intent(inout) :: nI(nel)
        integer, intent(in) :: nopen
        integer, intent(out) :: ncsf
        integer, intent(in) :: S
        logical, intent(in) :: tForceChange

        integer :: yamas (0:get_num_csfs(nopen, S), nopen), num
        real(dp) :: r

        if (nopen == 0) then
            ncsf = 0
            return
        endif

        ! Generate the Yamanouchi Symbols
        ncsf = size(yamas(:,1))-1
        call csf_get_yamas (nopen, S, yamas(1:,:), ncsf)

        if (tForceChange .and. ncsf > 1) then
            call get_csf_yama (nI, yamas(0,:), nopen)
        endif

        ! Pick and apply a random one
        do while (.true.)
            r = genrand_real2_dSFMT()
            num = int(r*ncsf) + 1
            if ((.not.tForceChange) .or. (ncsf<2) .or. &
                any(yamas(num,:) /= yamas(0,:))) then

                call csf_apply_yama (nI, yamas(num, :))
                exit
            endif
        enddo
    end subroutine

    subroutine csf_apply_yama (NI, csf)

        ! Apply a Yamanouchi symbol to a csf
        !
        ! In:    csf - The Yamanouchi symbol to apply
        ! InOut: nI  - The CSF to apply a Yamanouchi symbol to.

        integer, intent(in), dimension(:) :: csf
        integer, intent(inout) :: NI(nel)
        integer i

        nI = iand(nI, csf_orbital_mask)
        nI = ibset(nI, csf_test_bit)
        do i=1,size(csf)
            if (csf(size(csf)-i+1) .eq. 1) then
                NI(nel-i+1) = ibset(NI(nel-i+1), csf_yama_bit)
            else
                NI(nel-i+1) = ibclr(NI(nel-i+1), csf_yama_bit)
            endif
        enddo
    end subroutine

    subroutine csf_apply_ms (NI, Ms, nopen)

        ! Apply a specified Ms value to the csf, where we apply Ms/2
        !
        ! 

        integer, intent(inout) :: NI(nel)
        integer, intent(in) :: nopen
        integer, intent(in) :: Ms
        integer i, ndown
        
        ndown = (nopen - MS)/2
        do i=1,ndown
            NI(nel-i+1) = ibclr(NI(nel-i+1), csf_ms_bit)
        enddo
        do i=ndown+1,nopen
            NI(nel-i+1) = ibset(NI(nel-i+1), csf_ms_bit)
        enddo
    end subroutine

   integer function csf_spin (csf)

        ! Calculates the total spin from a csf
        !
        ! In:  csf - The Yamanouchi symbol to consider.
        ! Ret:     - 2 * the total spin

        integer, intent(in), dimension(:) :: csf
        integer i

        csf_spin = 0
        do i=1,size(csf)
            csf_spin = csf_spin + (1 - 2*(csf(i)-1))
        enddo
    end function

    integer function num_S (nopen)

        ! Calculates the number of possible S values (not their degeneracy)
        ! given a certain number of open shell electrons
        !
        ! In: nopen - Number of open shell electrons

        integer, intent(in) :: nopen
        num_S = (nopen+2)/2
    end function

    real(dp) pure function csf_coeff (csf, dorder, nopen)

        ! Calculate the coefficients for each determinant contained in the
        ! CSF. These are calculated as the product of Clebsch-Gordon coeffs.
        ! working through the tree electron-by-electron. Each coeff. depends
        ! on the current total spin in the csf, the current total spin in
        ! the determinant and the spin of the current e-/posn being considered
        ! in either the determinant or the csf.
        !
        ! In:  csf    - The Yamanouchi symbol to consider
        !      dorder - The list of alpha/beta for each spin orbital it det.
        !      nopen  - Number of unpaired electrons.
        ! Ret:        - Coefficient of the determinant represented by dorder,
        !               in the CSF represented by csf.

        integer, intent(in), dimension(:) :: csf, dorder
        integer, intent(in) :: nopen

        real(dp) :: S, M, scur, mcur, clb
        integer :: i

        S=0
        M=0
        csf_coeff = 1
        do i=1,nopen
            scur = -(csf(i)-1.5)
            mcur = -(real(dorder(i))-0.5)
            S = S + scur
            M = M + mcur

            clb = clbgrdn(S, M, scur, mcur)
            csf_coeff = csf_coeff * clb
            if (clb == 0) exit            
        enddo
    end function

    ! The Clebsch-Gordon coefficient with total spin S, projected
    ! component of spin M, with state with spin change 
    ! spin=(+/-)0.5 and projected spin change sig=(+/-)0.5
    real(dp) pure function clbgrdn(S,M,spin,sig)
        real(dp), intent(in) :: S, M, spin, sig
        !if (abs(spin).ne.0.5) then
        !    call stop_all ("clbgrdn","Spin incorrect")
        !endif
        if (spin.gt.0) then
            clbgrdn = sqrt((S+2*sig*M)/(2*S))
        else
            clbgrdn = -2*sig*sqrt((S+1-2*sig*M)/(2*(S+1)))
        endif
    end function

    subroutine write_yama (nunit, yama, lTerm)

        ! Write a Yamanouchi symbol to the output specified by nunit.
        ! if lTerm==.true. then terminate the line at the end.

        integer, intent(in) :: nunit
        integer, intent(in), dimension(:) :: yama
        logical, intent(in) :: lTerm
        integer i

        write (nunit,'("(")',advance='no')
        do i=1,size(yama)
            write(nunit,'(i1)',advance='no') yama(i)
        enddo
        write(nunit,'(")")',advance='no')
        if (lTerm) write(nunit,*)
    end subroutine


    function determine_s_converged (iluts, coeffs) result(s_final)

        ! Determine the S value for a converged set of determinants.
        !
        ! --> There are no determinants in a final wavefunction with a number
        !     of unpaired electrons lower than 2S (as this wavefunction needs
        !     to be equivalent by rotation to a wavefunction with S = Ms).
        !
        ! In:  iluts   - An array of bit representations of determinants to
        !                use.
        !      coeffs  - (optional) The coefficients to use for the 
        !                determinants provided. If not provided, use the sign
        !                value from the iluts.
        ! Ret: s_final - 2*S (i.e. return an integer value).

        integer(n_int), intent(in) :: iluts(:,:)
        real(dp), intent(in), optional :: coeffs(:)
        integer :: s_final

        integer :: i, sgn(lenof_sign), nopen, nopen_min
        real(dp) :: threshold, c


        ! If we are using sign values, then consider the cutoff value to be
        ! equivalent to the initiator threshold
        if (present(coeffs)) then
            threshold = real(InitiatorWalkNo, dp)
        else
            threshold = 1e-5
        endif

        ! We want to find the minimum number of unpaired electrons
        nopen_min = nel

        ! Loop over all the provided dets.
        do i = lbound(iluts, 2), ubound(iluts, 2)

            nopen = count_open_orbs(iluts(:,i))

            if (nopen < nopen_min) then

                if (present(coeffs)) then
                    c = coeffs(i)
                else
                    call extract_sign (iluts(:,i), sgn)
                    c = ARR_ABS(sgn)
                endif

                if (c > threshold) nopen_min = nopen

            endif
        enddo

        ! S = (nopen_min / 2). We multiply this value by 2 to represent it
        ! within an integer
        !
        ! --> return nopen_min

        s_final = nopen_min

    end function

    subroutine extract_dorder (nI, dorder, nopen)

        ! Return the dorder and the number of unpaired electrons in a
        ! determinant
        !
        ! In:  nI     - integer, ordered representation of det
        ! Out: dorder - Ordered list of alpha/beta for unpaired electrons.
        !               alpha == 0, beta == 1
        !      nopen  - Number of unpaired electrons in nI

        integer, intent(in) :: nI(nel)
        integer, intent(out) :: dorder(nel), nopen

        integer :: i

        nopen = 0
        i = 1
        do while (i <= nel)
            if (is_beta(nI(i)) .and. i < nel) then
                if (is_in_pair(nI(i), nI(i+1))) then
                    i = i + 2
                    cycle
                endif
            endif

            nopen = nopen + 1
            if (is_alpha(nI(i))) then
                dorder(nopen) = 0
            else
                dorder(nopen) = 1
            endif

            i = i + 1
        enddo

    end subroutine



end module

