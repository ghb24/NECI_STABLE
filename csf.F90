#include "macros.h"

! A new implementation file for csfs
module csf
    use systemdata, only: nel, brr, ecore, alat, nmsh, nbasismax, G1, nbasis
    use systemdata, only: NIfY, LMS, NIfTot, NIfD
    use memorymanager, only: LogMemAlloc, LogMemDealloc
    use integralsdata, only: umat, fck, nmax
    use HElem
    use mt95, only: genrand_real2
    use sltcnd_csf_mod, only: sltcnd_csf
    use DetBitOps, only: EncodeBitDet, FindBitExcitLevel
    use DetBitOps, only: get_bit_excitmat_op_ind
    use csf_data
    use, intrinsic :: ieee_arithmetic

    implicit none

    ! Non-modularised functions (sigh)
    interface
        real*8 pure function choose(N,R)
            integer, intent(in) :: N,R
        end function
        logical function int_arr_eq (a, b, len)
            integer, intent(in), dimension(:) :: a, b
            integer, intent(in), optional :: len
        end function
        function GetHElement3_wrapper(NI,NJ,iC)
            use HElem
            use Systemdata, only: nEl
            INTEGER NI(nEl),NJ(nEl),iC
            type(HElement) GetHElement3_wrapper
        end function
        subroutine writedet(nunit,ni,nel,lterm)
            integer nunit,nel,ni(nel)
            logical lterm
        end subroutine
    end interface
    
contains
    ! The number of determinants given a specified number of csfs (note that
    ! if nopen == 0, ie ncsf==0, there is still a closed shell determinant
    ! legitimately exists.
    integer function num_csf_dets (ncsf)
        integer, intent(in) :: ncsf
        if (ncsf == 0) then
            num_csf_dets = 1
        else
            num_csf_dets = ncsf
        endif
    end function

    function CSFGetHelement (nI, nJ) result (hel_ret)
        integer, intent(in) :: nI(nel), nJ(nel)
        integer :: ilutI(0:NIfTot), iLutJ(0:NIfTot)
        type(HElement) :: hel_ret

        call EncodeBitDet (nI, iLutI)
        call EncodeBitDet (nJ, iLutJ)

        hel_ret = CSFGetHelement_bit(nI, nJ, iLutI, iLutJ)
    end function

    ! This is (intentionally) a BRUTE FORCE way of calculating this.
    ! We would like to see if there is a nicer way of doing this using
    ! the representation matrices of the permutations which we are able 
    ! to calculate.
    ! TODO: Can probably do a lot of this stuff faster with bit dets...
    function CSFGetHelement_bit(NI, NJ, iLutI, iLutJ) result(hel_ret)
        implicit none
        integer, intent(in) :: NI(nel), NJ(nel)
        integer nopen(2), nclosed(2), nup(2), ndets(2), i, j, det
        real*8 S(2), Ms(2)
        type(HElement) Hel, sum1, hel_ret

        integer, intent(in) :: iLutI(0:NIfTot), iLutJ(0:NIfTot)
        integer, dimension (:), allocatable :: yama1, yama2
        real*8, dimension (:), allocatable :: coeffs1, coeffs2
        integer, dimension (:,:), allocatable :: dets1, dets2
        !integer, dimension(:,:), allocatable :: dets1a, dets2a
        integer, dimension(:,:), allocatable :: ilut1, ilut2
        integer tagYamas(2)/0,0/, tagCoeffs(2)/0,0/, tagDets(2)/0,0/
        integer tagILuts(2)/0,0/
        integer :: IC

        character(*), parameter :: this_routine = 'CSFGetHelement'

        ! If the determinants differ by more than 2 spacial orbitals, then 
        ! they cannot contribute.
        IC = FindBitExcitLevel (iLutI, iLutJ)
        if (IC > 2) then
            hel_ret = HElement(0)
            return
        endif

        call get_csf_data(NI, nel, nopen(1), nclosed(1), S(1), Ms(1))
        call get_csf_data(NJ, nel, nopen(2), nclosed(2), S(2), Ms(2))

        ! If S or Ms are not consistent, then return 0
        if ((S(1).ne.S(2))) then ! .or. (Ms(1).ne.Ms(2))) then
            hel_ret = HElement(0) 
            return
        endif
        Ms(1) = minval(S)
        Ms(2) = Ms(1)

        ! Get electronic details
        ! Using S instead of Ms to calculate nup, as this has the fewest
        ! determinants, and the Ms=S case is degenerate.
        nup(1) = (nopen(1) + 2*S(1))/2
        ndets(1) = int(choose(nopen(1),nup(1)))
        nup(2) = (nopen(2) + 2*S(2))/2
        ndets(2) = int(choose(nopen(2),nup(2)))

        ! Allocate as required
        ! TODO: catch ierr if it fails.
        allocate (coeffs1(ndets(1)), coeffs2(ndets(2)))
        allocate (dets1(ndets(1), nel), dets2(ndets(2), nel), &
                   !dets1a(ndets(1), nel), dets2a(ndets(2), nel), &
                  ilut1(ndets(1),0:NIfTot), ilut2(ndets(2),0:NIfTot))
        if (nopen(1) > 0) then
            allocate(yama1(nopen(1)))
            call LogMemAlloc ('yama1',size(yama1),4,this_routine,tagYamas(1))
        endif
        if (nopen(2) > 0) then
            allocate(yama2(nopen(2)))
            call LogMemAlloc ('yama2',size(yama2),4,this_routine,tagYamas(2))
        endif
        call LogMemAlloc ('coeffs1',size(coeffs1),8,this_routine,tagCoeffs(1))
        call LogMemAlloc ('coeffs2',size(coeffs2),8,this_routine,tagCoeffs(2))
        call LogMemAlloc ('dets1',size(dets1),4,this_routine,tagDets(1))
        call LogMemAlloc ('dets2',size(dets2),4,this_routine,tagDets(2))
        call LogMemAlloc ('ilut1',size(ilut1),4,this_routine,tagIluts(1))
        call LogMemAlloc ('ilut2',size(ilut2),4,this_routine,tagIluts(2))
        
        ! Calculate all possible permutations to construct determinants
        ! (Where 0=alpha, 1=beta when generating NI, NJ below)
        call csf_get_dets(nopen(1), nup(1), ndets(1), nel, dets1)
        if ((nopen(1).eq.nopen(2)) .and. (nup(1).eq.nup(2))) then
            dets2 = dets1
        else
            call csf_get_dets(nopen(2), nup(2), ndets(2), nel, dets2)
        endif

        ! Extract the Yamanouchi symbols from the CSFs
        call get_csf_yama (NI, yama1)
        call get_csf_yama (NJ, yama2)

        ! Get the coefficients
        do det=1,ndets(1)
            coeffs1(det) = csf_coeff(yama1,dets1(det,nclosed(1)+1:nel),&
                                     nopen(1))
        enddo
        do det=1,ndets(2)
            coeffs2(det) = csf_coeff(yama2,dets2(det,nclosed(2)+1:nel),&
                                     nopen(2))
        enddo

        ! Generate determinants from spatial orbitals specified in NI, NJ
        do det = 1,ndets(1)
            dets1(det,1:nclosed(1)) = iand(NI(1:nclosed(1)), csf_orbital_mask)
            dets1(det,nclosed(1)+1:nel) = &
                    csf_alpha_beta(NI(nclosed(1)+1:nel), &
                                   dets1(det,nclosed(1)+1:nel))
            call EncodeBitDet (dets1(det,:), ilut1(det,:))
        enddo
        do det = 1,ndets(2)
            dets2(det,1:nclosed(2)) = iand(NJ(1:nclosed(2)), csf_orbital_mask)
            dets2(det,nclosed(2)+1:nel) = &
                    csf_alpha_beta(NJ(nclosed(2)+1:nel),&
                                   dets2(det,nclosed(2)+1:nel))
            call EncodeBitDet (dets2(det,:), ilut2(det,:))
        enddo

        hel_ret = HElement(0)
        do i=1,ndets(1)
            if (coeffs1(i) /= 0) then
                sum1 = Helement(0)
                do j=1,ndets(2)
                    if (coeffs2(j) /= 0) then
                        Hel = sltcnd_csf (dets1(i,:), dets2(j,:), &
                                          ilut1(i,:), ilut2(j,:))
                        sum1 = sum1 + Hel * HElement(coeffs2(j))
                    endif
                enddo
                hel_ret = hel_ret + sum1*HElement(coeffs1(i))
            endif
        enddo

        ! Deallocate for cleanup
        deallocate (coeffs1, coeffs2, dets1, dets2, ilut1, ilut2)
        if (allocated(yama1)) then
            deallocate(yama1)
            call LogMemDealloc (this_routine, tagYamas(1))
        endif
        if (allocated(yama2)) then
            deallocate(yama2)
            call LogMemDealloc (this_routine, tagYamas(2))
        endif
        call LogMemDealloc (this_routine, tagCoeffs(1))
        call LogMemDealloc (this_routine, tagCoeffs(2))
        call LogMemDealloc (this_routine, tagDets(1))
        call LogMemDealloc (this_routine, tagDets(2))
        call LogMemDealloc (this_routine, tagILuts(1))
        call LogMemDealloc (this_routine, tagILuts(2))
    end function

    function CSFGetHelement_faster(nI, nJ) result(hel_ret)
        
        ! Calculate the H-matrix element between two CSFs (nI, nJ)
        !
        ! In:  nI, nJ   - The determinants to consider
        ! Ret: hel_ret  - The H-matrix element

        integer, intent(in) :: NI(nel), NJ(nel)
        type(HElement) :: hel_ret

        integer :: nopen(2), nclosed(2), nup(2), ndets(2)
        integer :: iLutI(0:NIfTot), iLutJ(0:NIfTot), IC
        ! Convert these to using integers
        real*8  :: S(2), Ms(2)

        character(*), parameter :: this_routine = 'CSFGetHelement_faster'


        ! TODO: We should be able to pass these through
        ! If the determinants differ by more than 2 spacial orbitals, then 
        ! they cannot contribute.
        call EncodeBitDet (nI, iLutI)
        call EncodeBitDet (nJ, iLutJ)

        IC = FindBitExcitLevel (iLutI, iLutJ)
        if (IC > 2) then
            hel_ret = HElement(0)
            return
        endif

        call get_csf_data(NI, nel, nopen(1), nclosed(1), S(1), Ms(1))
        call get_csf_data(NJ, nel, nopen(2), nclosed(2), S(2), Ms(2))

        ! If S or Ms are not consistent, then return 0
        if ((S(1).ne.S(2))) then ! .or. (Ms(1).ne.Ms(2))) then
            hel_ret = HElement(0) 
            return
        endif

        ! Use the maximal Ms value that we can (fewest determinants required)
        Ms(1) = minval(S)
        Ms(2) = Ms(1)

        ! Get electronic details
        ! Using S instead of Ms to calculate nup, as this has the fewest
        ! determinants, and the Ms=S case is degenerate.
        nup(1) = (nopen(1) + 2*S(1))/2
        ndets(1) = int(choose(nopen(1),nup(1)))
        nup(2) = (nopen(2) + 2*S(2))/2
        ndets(2) = int(choose(nopen(2),nup(2)))

        hel_ret = get_csf_helement_local (nI, nJ, iLutI, iLutJ, nopen, &
                                          nclosed, nup, ndets, IC)
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
        integer, intent(in) :: iLutI(0:NIfTot), iLutJ(0:NIfTot)
        integer, intent(in) :: nup(2), ndets(2), IC
        type(HElement) :: hel_ret

        integer :: yama1(nopen(1)), yama2(nopen(2))
        real*8 :: coeffs1(ndets(1)), coeffs2(ndets(2))
        integer :: dets1(nel, ndets(1)), dets2(nel, ndets(2))
        integer :: ilut1(0:NIfTot,ndets(1)), ilut2(0:NIfTot,ndets(2))

        integer det, i, j
        type(HElement) :: sum1, Hel

        if (IC == 2) then
            hel_ret = get_csf_helement_2 (nI, nJ, iLutI, iLutJ, nopen,  &
                                    nclosed, nup, ndets, dets1, dets2,&
                                    yama1, yama2, coeffs1, coeffs2, ilut1, &
                                    ilut2)
            return
        endif

        ! TODO: IC==0 symmetry
        ! If IC==0, then only differ by yama --> dets are the same. Have some
        ! symmetry on the summation --> half the number of terms!

        ! Calculate all possible permutations to construct determinants
        ! (Where 0=alpha, 1=beta when generating NI, NJ below)
        call csf_get_dets_reverse (nopen(1), nup(1), ndets(1), nel, dets1)
        if ((nopen(1).eq.nopen(2)) .and. (nup(1).eq.nup(2))) then
            dets2 = dets1
        else
            call csf_get_dets_reverse (nopen(2), nup(2), ndets(2), nel, dets2)
        endif

        ! Extract the Yamanouchi symbols from the CSFs
        call get_csf_yama (NI, yama1)
        call get_csf_yama (NJ, yama2)

        ! Get the coefficients
        do det=1,ndets(1)
            coeffs1(det) = csf_coeff(yama1,dets1(nclosed(1)+1:nel,det),&
                                     nopen(1))
        enddo
        do det=1,ndets(2)
            coeffs2(det) = csf_coeff(yama2,dets2(nclosed(2)+1:nel,det),&
                                     nopen(2))
        enddo

        ! Generate determinants from spatial orbitals specified in NI, NJ
        do det = 1,ndets(1)
            dets1(1:nclosed(1),det) = iand(NI(1:nclosed(1)), csf_orbital_mask)
            dets1(nclosed(1)+1:nel,det) = &
                    csf_alpha_beta(NI(nclosed(1)+1:nel), &
                                   dets1(nclosed(1)+1:nel,det))
            call EncodeBitDet (dets1(:,det), ilut1(:,det))
        enddo
        do det = 1,ndets(2)
            dets2(1:nclosed(2),det) = iand(NJ(1:nclosed(2)), csf_orbital_mask)
            dets2(nclosed(2)+1:nel,det) = &
                    csf_alpha_beta(NJ(nclosed(2)+1:nel),&
                                   dets2(nclosed(2)+1:nel,det))
            call EncodeBitDet (dets2(:,det), ilut2(:,det))
        enddo

        ! TODO: implement symmetry if NI,NJ are the same except for yama
        hel_ret = HElement(0)
        do i=1,ndets(1)
            if (coeffs1(i) /= 0) then
                sum1 = Helement(0)
                do j=1,ndets(2)
                    if (coeffs2(j) /= 0) then
                        Hel = sltcnd_csf (dets1(:,i), dets2(:,j), &
                                          ilut1(:,i), ilut2(:,j))
                        sum1 = sum1 + Hel * HElement(coeffs2(j))
                    endif
                enddo
                hel_ret = hel_ret + sum1*HElement(coeffs1(i))
            endif
        enddo
    end function

    function get_csf_helement_2 (nI, nJ, iLutI, iLutJ, nopen, nclosed, &
                                 nup, ndets, dets1, dets2, yama1, yama2, &
                                 coeffs1, coeffs2, ilut1, ilut2) &
                                 result(hel_ret)

        ! Given a case where IC == 2, calculate the helement. Thus we can make
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

        integer, intent(in) :: nI(nel), nJ(nel), nopen(2), nclosed(2)
        integer, intent(in) :: iLutI(0:NIfTot), iLutJ(0:NIfTot)
        integer, intent(in) :: nup(2), ndets(2)
        ! TODO: convert these to integers (ie 2S, 2Ms)
        integer, intent(out) :: yama1(nopen(1)), yama2(nopen(2))
        real*8, intent(out) :: coeffs1(ndets(1)), coeffs2(ndets(2))
        integer, intent(inout) :: dets1(nel, ndets(1)), dets2(nel, ndets(2))
        integer, intent(out) :: ilut1(0:NIfTot,ndets(1)), &
                                ilut2(0:NIfTot,ndets(2))
        type(HElement) :: hel_ret, hel, sum1

        integer :: nsing_delta(2), det, i, j, k
        integer :: ex_id(4,2), ms1(ndets(1)), ms2(ndets(2))
        logical :: dets_change1(ndets(1)), dets_change2(ndets(2))

        ! count the number of singles which differ between the two dets, and
        ! get their indices in the open section.
        call get_bit_open_unique_ind (iLutI, iLutJ, ex_id, nsing_delta, 2)

        ! Calculate all possible permutations to construct determinants
        ! (Where 0=alpha, 1=beta when generating NI, NJ below)
        call flush(6)
        call csf_get_dets_ind (nopen(1), nup(1), ndets(1), nel, &
                               nsing_delta(1), ex_id(:,1), dets1, ms1)
        call csf_get_dets_ind (nopen(2), nup(2), ndets(2), nel, &
                               nsing_delta(2), ex_id(:,2), dets2, ms2)

        ! TODO: can we do this in the previous bit?
        ! Mark each permutation for which each of the differing orbitals are
        ! permuted in a canonical order --> the end of the set which mixes.
        call mark_change_2 (dets_change1, dets1, nclosed(1), ndets(1), &
                            nsing_delta(1), ex_id(:,1))
        call mark_change_2 (dets_change2, dets2, nclosed(2), ndets(2), &
                            nsing_delta(2), ex_id(:,2))

        ! Extract the Yamanouchi symbols from the CSFs
        call get_csf_yama (nI, yama1)
        call get_csf_yama (nJ, yama2)

        ! Get the coefficients
        do det=1,ndets(1)
            coeffs1(det) = csf_coeff(yama1,dets1(nclosed(1)+1:nel,det),&
                                     nopen(1))
        enddo
        do det=1,ndets(2)
            coeffs2(det) = csf_coeff(yama2,dets2(nclosed(2)+1:nel,det),&
                                     nopen(2))
        enddo

        ! Generate determinants from spatial orbitals specified in NI, NJ
        do det = 1,ndets(1)
            if (coeffs1(det) /= 0) then
                dets1(1:nclosed(1),det) = nI(1:nclosed(1))
                dets1(nclosed(1)+1:nel,det) = &
                        csf_alpha_beta(nI(nclosed(1)+1:nel), &
                                       dets1(nclosed(1)+1:nel,det))
                call EncodeBitDet (dets1(:,det), ilut1(:,det))
            endif
        enddo
        do det = 1,ndets(2)
            if (coeffs2(det) /= 0) then
                dets2(1:nclosed(2),det) = nJ(1:nclosed(2))
                dets2(nclosed(2)+1:nel,det) = &
                        csf_alpha_beta(nJ(nclosed(2)+1:nel),&
                                       dets2(nclosed(2)+1:nel,det))
                call EncodeBitDet (dets2(:,det), ilut2(:,det))
            endif
        enddo

        ! TODO: neat tricks using excitmat to avoid call to sltcnd
        j = 1
        i = 1
        hel_ret = HElement(0)
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
            
            if (coeffs1(i) /= 0) then
                sum1 = HElement(0)
                k = j
                do while (k <= ndets(2))
                    if (coeffs2(k) /= 0) then
                        ! Optimise this away.
                        Hel = sltcnd_csf (dets1(:,i), dets2(:,k), &
                                          ilut1(:,i), ilut2(:,k))
                        sum1 = sum1 + Hel*HElement(coeffs2(k))
                    endif

                    k = k + 1
                    if (dets_change2(k-1)) exit
                enddo
                hel_ret = hel_ret + sum1*HElement(coeffs1(i))
            endif

            if (dets_change1(i)) then
                do while (j <= ndets(2))
                    j = j + 1
                    if (dets_change2(j-1)) exit
                enddo
            endif

            i = i + 1
        enddo

    end function

!    function get_csf_helement_local_bit (nI, nJ, iLutI, iLutJ, nopen, nclosed, S,&
!                                     Ms, nup, ndets) result(hel_ret)
!        integer, intent(in) :: nI(nel), nJ(nel), nopen(2), nclosed(2)
!        integer, intent(in) :: iLutI(0:NIfTot), iLutJ(0:NIfTot)
!        integer, intent(in) :: nup(2), ndets(2)
!        ! TODO: convert these to integers (ie 2S, 2Ms)
!        real*8, intent(in) :: S(2), Ms(2)
!        type(HElement) :: hel_ret
!
!        integer :: yama1(NIfY), yama2(NIfY)
!        real*8 :: coeffs1(ndets(1)), coeffs2(ndets(2))
!        integer :: perm1(NIfY, ndets(1)), perm2(NIfY, ndets(2))
!
!        integer det, i, j
!
!        ! Calculate all possible permutations to construct determinants
!        ! (Where 0=alpha, 1=beta when generating NI, NJ below)
!        call csf_get_bit_perm (nopen(1), nup(1), ndets(1), nel, perm1)
!        if ((nopen(1).eq.nopen(2)) .and. (nup(1).eq.nup(2))) then
!            perm2 = perm1
!        else
!            call csf_get_bit_perm(nopen(2), nup(2), ndets(2), nel, perm2)
!        endif
!
!        ! Extract the Yamanouchi symbols from the CSFs
!        yama1 = iLutI (NIfD+1:NIfD+NIfY)
!        yama2 = iLutJ (NIfD+1:NIfD+NIfY)
!
!        ! Get the coefficients
!        do det=1,ndets(1)
!            coeffs1(det) = csf_bit_coeff(yama1, perm1(:,det), nopen(1))
!        enddo
!        do det=1,ndets(2)
!            coeffs2(det) = csf_bit_coeff(yama2, perm2(:,det), nopen(2))
!        enddo
!
!
!
!        ! TODO: implement symmetry if NI,NJ are the same except for yama
!        !hel_ret = HElement(0)
!        !do i=1,ndets(1)
!        !    if (coeffs1(i) /= 0) then
!        !        sum1 = Helement(0)
!        !        do j=1,ndets(2)
!        !            if (coeffs2(j) /= 0) then
!        !                Hel = sltcnd_csf (dets1(i,:), dets2(j,:), &
!        !                                  ilut1(i,:), ilut2(j,:))
!        !                sum1 = sum1 + Hel * HElement(coeffs2(j))
!        !            endif
!        !        enddo
!        !        hel_ret = hel_ret + sum1*HElement(coeffs1(i))
!        !    endif
!        !enddo
!
!    end function
    
    ! TODO: comment
    subroutine mark_change_2 (dets_change, dets, nclosed, ndets, nopen_sing,&
                              ex_id)
        integer, intent(in) :: nclosed, ndets, nopen_sing
        integer, intent(in) :: dets(nel, ndets), ex_id(4)
        logical, intent(out) :: dets_change(ndets)
        logical :: bChange
        integer :: i, j

        dets_change = .true.
        if (nopen_sing == 0) return

        do i=1,ndets
            bChange = .false.
            do j=1,nopen_sing
                if (.not. bChange .and. dets(ex_id(j)+nclosed,i) == 0) &
                    bChange = .true.

                if (bChange .and. dets(ex_id(j)+nclosed,i) /= 0) then
                    dets_change(i) = .false.
                    exit
                endif
            enddo                    
        enddo

    end subroutine

    subroutine csf_get_bit_perm (nopen, nup, ndets, ilut)
    
        ! TODO: comment and reverse ordering as below.
        
        integer, intent(in) :: ndets, nup, nopen
        integer, intent(out) :: ilut(NIfY,ndets)
        integer :: i, det, comb(nup)
        logical :: bInc

        if (nopen == 0) return

        forall (i=1:nup) comb(i) = i
        ilut = 1
        do det = 1, ndets
            forall (i=1:nup) ilut((comb(i)-1)/32,det) = &
                            ibclr(ilut((comb(i)-1)/32,det), mod(comb(i)-1,32))

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

    subroutine csf_get_dets_reverse (nopen, nup, ndets, nel, dets)

        ! Fill the last nopen electrons of each determinant with 0 (alpha) or
        ! 1 (beta) in all possible permutations with nup alpha electrons.

        integer, intent(in) :: ndets, nup, nopen, nel
        integer, intent(out) :: dets (nel,ndets)
        integer comb(nup), i, j
        logical bInc

        if (nopen.eq.0) return

        forall (i=1:nup) comb(i) = nopen-i+1
        ! TODO: fix ordering
        dets(nel-nopen+1:,:) = 1
        do i=1,ndets
            forall (j=1:nup) dets(nel-nopen+comb(j),i) = 0
            do j=1,nup
                bInc = .false.
                if (j == nup) then
                    bInc = .true.
                else if (j < nup) then
                    if (comb(j+1) /= comb(j) - 1) bInc = .true.
                endif

                if (bInc) then
                    comb(j) = comb(j) - 1
                    exit
                else
                    comb(j) = nopen - j + 1
                endif
            enddo
        enddo
    end subroutine

    subroutine csf_get_dets (nopen, nup, ndets, nel, dets)

        ! Fill the last nopen electrons of each determinant with 0 (alpha) or
        ! 1 (beta) in all possible permutations with nup alpha electrons.

        integer, intent(in) :: ndets, nup, nopen, nel
        integer, intent(out) :: dets (ndets,nel)
        integer comb(nup), i, j
        logical bInc

        if (nopen.eq.0) return

        forall (i=1:nup) comb(i) = i !nopen-i+1
        ! TODO: fix ordering
        dets(:,nel-nopen+1:) = 1
        do i=1,ndets
            forall (j=1:nup) dets(i,nel-nopen+comb(j)) = 0
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
                    comb(j) = j!nopen-j+1
                endif
            enddo
        enddo
    end subroutine

    subroutine csf_get_dets_ind (nopen, nup, ndets, nel, nuniq, ex_id, dets, &
                                 ms)

        ! TODO: comment
        ! Fill the last nopen electrons of each determinant with 0 (alpha) or
        ! 1 (beta) in all possible permutations with nup alpha electrons.

        integer, intent(in) :: ndets, nup, nopen, nel, nuniq, ex_id(*)
        integer, intent(out) :: dets (nel,ndets), ms(ndets)
        integer comb(nup), i, j, id(nopen), pos, posu
        logical bInc

        if (nopen.eq.0) then
            ms = 0
            !print*, ' 0|'
            return
        endif

        ! Generate mapping of positions to cause the unique indices to vary as
        ! fast as possible
        id(1:nuniq) = ex_id(1:nuniq)
        pos = nuniq+1
        posu = 1
        do i=1,nopen
            if (posu > nuniq) exit
            if (i == ex_id(posu)) then
                posu = posu + 1
            else
                id(pos) = i
                pos = pos + 1
            endif
        enddo
        if (i <= nopen) forall (j=i:nopen) id(j) = j

        ! Calculate permutations and place in dets
        forall (i=1:nup) comb(i) = i
        dets(nel-nopen+1:,:) = 1
        do i=1,ndets
            forall (j=1:nup) dets(nel-nopen+id(comb(j)),i) = 0
            ms(i) = nup - sum(dets(nel-nopen+id(1:nuniq),i))
            ! TODO: remove
            !write (6,'(i3,"|")', advance='no') ms(i)
            !print*, dets(nel-nopen+id(1:nuniq),i)
            !print*, dets(nel-nopen+id(nuniq+1:nopen),i)
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

    ! For all of the open shell electrons, generate a list where
    ! 0=alpha, 1=beta for the specified determinant.
    function csf_get_dorder (NI, nel, nopen)
        integer csf_get_dorder(nopen)
        integer, intent(in) :: NI(nel)
        integer, intent(in) :: nopen, nel
        integer nclosed

        nclosed = nel - nopen
        csf_get_dorder = -(mod(NI(nclosed+1:nel),2)-1)
    end function


    subroutine csf_get_yamas (nopen, sfinal, yama, ncsf_max)
        real*8, intent(in) :: sfinal
        integer, intent(in) :: nopen, ncsf_max
        integer, intent(out) :: yama (ncsf_max, nopen)
        real*8 spin (ncsf_max, nopen)
        integer npos, csf, ncsf, ncsf_next

        if (nopen == 0) return

        spin(1,nopen) = sfinal
        ncsf = 1
        ncsf_next = ncsf
        do npos = nopen, 2, -1
            do csf=1,ncsf
                if (2*spin(csf,npos) .lt. npos) then
                    spin(csf,npos-1) = spin(csf,npos) + 0.5
                    yama(csf,npos) = 2
                    if (spin(csf,npos) .ne. 0) then
                        ncsf_next = ncsf_next + 1
                        if (ncsf_next .gt. ncsf_max) exit
                        !spin(ncsf_next,npos:nopen) = spin(csf,npos:nopen)
                        yama(ncsf_next,npos+1:nopen) = yama(csf,npos+1:nopen)
                        spin(ncsf_next,npos-1) = spin(csf,npos) - 0.5
                        yama(ncsf_next,npos) = 1
                    endif
                else
                    spin(csf,npos-1) = spin(csf,npos) - 0.5
                    yama(csf,npos) = 1
                endif
            enddo
            ncsf = ncsf_next
            if (ncsf .gt. ncsf_max) exit
        enddo
        yama(:,1) = 1
    end subroutine
        
    ! Convert num (a member of CI) to an alpha or beta spin orbital where
    ! det==0 --> alpha
    integer elemental function csf_alpha_beta (num, det)
        integer, intent(in) :: num, det

        csf_alpha_beta = iand(num, csf_orbital_mask) - 1
        if (det .eq. 1) then
            csf_alpha_beta = ior(csf_alpha_beta, 1)
        else
            csf_alpha_beta = iand(csf_alpha_beta,Z'fffffffe')
        endif
        csf_alpha_beta = csf_alpha_beta + 1
    end function

    ! Extract the Yamanouchi symbol from the supplied CSF
    ! This assumes that the passed yama array is the correct size, if
    ! it is too small, then the symbol will be truncated.
    subroutine get_csf_yama(NI, yama)
        integer, intent(in), dimension(:) :: NI
        integer, intent(out), dimension(:) :: yama
        integer i, nopen
        logical open_shell

        nopen = 0
        open_shell = .false.
        do i=1,size(NI)
            if ( (.not.open_shell) .and. btest(NI(i), csf_yama_bit)) then
                open_shell = .true.
            endif

            if (open_shell) then
                nopen = nopen + 1
                if (btest(NI(i), csf_yama_bit)) then
                    yama(nopen) = 1
                else
                    yama(nopen) = 2
                endif
                if (nopen == size(yama)) exit
            endif
        enddo
    end subroutine

    subroutine get_csf_bit_yama (nI, yama)
        
        ! Extract the Yamanouchi symbol from the supplied CSF as part of
        ! a bit determinant - the same as part of the bit det obtained by
        ! EncodeBitDet.

        integer, intent(in) :: nI(nel)
        integer, intent(out) :: yama(NIfY)
        integer :: i, pos, bit, nopen

        nopen = 0
        yama = 0
        do i=1,nel
            ! Only consider open electrons. The first has csf_yama_bit set.
            if (nopen > 0 .or. btest(nI(i), csf_yama_bit)) then
                pos = 1 + nopen / 32
                bit = mod(nopen, 32)

                if (btest(nI(i), csf_yama_bit)) &
                    yama(pos) = ibset(yama(pos), bit)
                nopen = nopen + 1
            endif
        enddo
    end subroutine

    subroutine csf_to_old_csf (nI, nJ)
        
        ! Convert the CSF to the old representation (Alex's rep.) for testing
        ! purposes.
        
        integer, intent(in) :: nI(nel)
        integer, intent(out) :: nJ(nel)
        integer :: i, nopen, nup, orb
        logical :: open_shell
        include 'csf.inc'

        if (.not. iscsf(nI)) then
            nJ = nI
            return
        endif

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



    ! Obtains the number of open shell electrons, the total spin
    ! and the Ms value for the specified csf.
    subroutine get_csf_data(NI, nel, nopen, nclosed, S, Ms)
        integer, intent(in) :: NI(nel), nel
        integer, intent(out) :: nopen, nclosed
        real*8, intent(out) :: S, Ms
        integer i
        logical open_shell

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
                    S = S + 0.5
                else
                    S = S - 0.5
                endif
                if (btest(NI(i), csf_ms_bit)) then
                    Ms = Ms + 0.5
                else
                    Ms = Ms - 0.5
                endif
            endif
        enddo
        nopen = nel - nclosed
    end subroutine

    ! Calculates the total number of CSFs possible for a system with
    ! nOpen unpaired electrons, and a total spin of S.
    ! This is the same as the number of available Serber functions
    integer pure function get_num_csfs (nOpen, S)
        integer, intent(in) :: nOpen
        real*8, intent(in) :: S
        integer :: S2
        S2 = 2*S

        if ((nopen < 0) .or. (mod(nOpen+S2, 2) /= 0))then
            get_num_csfs = 0
        else
            get_num_csfs = (2*S2 + 2) * choose(nOpen, (nOpen+S2)/2)
            get_num_csfs = get_num_csfs / (nOpen + S2 + 2)
        endif
    end function

    ! TODO: This can be optimised (don't need to generate them all)
    ! TODO: Generate random by random branching perhaps (lots of genrands...)
    subroutine csf_apply_random_yama (nI, nopen, S, ncsf, tForceChange)
        integer, intent(inout) :: nI(nel)
        integer, intent(in) :: nopen
        integer, intent(out) :: ncsf
        real*8, intent(in) :: S
        logical, intent(in) :: tForceChange
        integer :: yamas (0:get_num_csfs(nopen, S), nopen), num
        real*8 :: r

        ! Generate the Yamanouchi Symbols
        ncsf = size(yamas(:,1))-1
        call csf_get_yamas (nopen, S, yamas(1:,:), ncsf)

        if (tForceChange .and. ncsf > 1) then
            call get_csf_yama (nI, yamas(0,:))
        endif

        ! Pick and apply a random one
        do while (.true.)
            call genrand_real2(r)
            num = int(r*ncsf) + 1
            if ((.not.tForceChange) .or. (ncsf<2) .or. &
                .not.int_arr_eq(yamas(num,:),yamas(0,:))) then

                call csf_apply_yama (nI, yamas(num, :))
                exit
            endif
        enddo
    end subroutine

    ! Apply a Yamanouchi symbol to a csf
    subroutine csf_apply_yama (NI, csf)
        integer, intent(in), dimension(:) :: csf
        integer, intent(inout) :: NI(nel)
        integer i

        NI = ibset(NI, csf_test_bit)
        do i=1,size(csf)
            if (csf(size(csf)-i+1) .eq. 1) then
                NI(nel-i+1) = ibset(NI(nel-i+1), csf_yama_bit)
            else
                NI(nel-i+1) = ibclr(NI(nel-i+1), csf_yama_bit)
            endif
        enddo
    end subroutine

    ! Apply a specified Ms value to the csf
    subroutine csf_apply_ms (NI, Ms, nopen)
        integer, intent(inout) :: NI(nel)
        integer, intent(in) :: nopen
        real*8, intent(in) :: Ms
        integer i, ndown
        
        ndown = (nopen - 2*MS)/2
        do i=1,ndown
            NI(nel-i+1) = ibclr(NI(nel-i+1), csf_ms_bit)
        enddo
        do i=ndown+1,nopen
            NI(nel-i+1) = ibset(NI(nel-i+1), csf_ms_bit)
        enddo
    end subroutine

    ! Calculates the total spin from a csf
    real*8 function csf_spin (csf)
        integer, intent(in), dimension(:) :: csf
        integer i

        csf_spin = 0
        do i=1,size(csf)
            csf_spin = csf_spin - (real(csf(i))-1.5)
        enddo
    end function

    ! Calculates the number of possible S values (not their degeneracy)
    ! given a certain number of open shell electrons
    integer function num_S (nopen)
        integer, intent(in) :: nopen
        num_S = (nopen+2)/2
    end function

    ! Calculate the coefficients for each determinant contained in the
    ! CSF. These are calculated as the product of Clebsch-Gordon coeffs.
    ! working through the tree electron-by-electron. Each coeff. depends
    ! on the current total spin in the csf, the current total spin in
    ! the determinant and the spin of the current e-/posn being considered
    ! in either the determinant or the csf.
    ! dorder = the ordered list of alpha/beta for each spin orbital in det.
    real*8 pure function csf_coeff (csf, dorder, nopen)
        integer, intent(in), dimension(:) :: csf, dorder
        integer, intent(in) :: nopen
        real*8 S, M, scur, mcur, clb
        integer i

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

    ! Calculate the coefficients for each determinant contained in the
    ! CSF. These are calculated as the product of Clebsch-Gordon coeffs.
    ! working through the tree electron-by-electron. Each coeff. depends
    ! on the current total spin in the csf, the current total spin in
    ! the determinant and the spin of the current e-/posn being considered
    ! in either the determinant or the csf.
    ! dorder = the ordered list of alpha/beta for each spin orbital in det.
    real*8 pure function csf_bit_coeff (csf, dorder, nopen)
        integer, intent(in) :: csf(NIfY), dorder(NIfY), nopen
        real*8 S, M, scur, mcur, clb
        integer i, bit, pos

        S=0
        M=0
        csf_bit_coeff = 1
        do i=0,nopen-1
            bit = mod(i,32)
            pos = i/32

            scur = -0.5
            mcur = 0.5
            if (ibset(csf(pos), bit)) scur = 0.5
            if (ibset(dorder(pos), bit)) scur = -0.5
            S = S + scur
            M = M + mcur

            clb = clbgrdn(S, M, scur, mcur)
            csf_bit_coeff = csf_bit_coeff * clb
            if (clb == 0) exit            
        enddo
    end function

    ! The Clebsch-Gordon coefficient with total spin S, projected
    ! component of spin M, with state with spin change 
    ! spin=(+/-)0.5 and projected spin change sig=(+/-)0.5
    real*8 pure function clbgrdn(S,M,spin,sig)
        real*8, intent(in) :: S,M,spin,sig
        !if (abs(spin).ne.0.5) then
        !    call stop_all ("clbgrdn","Spin incorrect")
        !endif
        if (spin.gt.0) then
            clbgrdn = sqrt((S+2*sig*M)/(2*S))
        else
            clbgrdn = -2*sig*sqrt((S+1-2*sig*M)/(2*(S+1)))
        endif
    end function

    ! Write a Yamanouchi symbol to output specified by nunit.
    ! if lTerm==.true. then terminate the line at the end.
    subroutine write_yama (nunit, yama, lTerm)
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
end module

