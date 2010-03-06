#include "macros.h"

module sltcnd_mod
    use SystemData, only: nel, nBasisMax, tExch, FCOUL, NIfTot, G1, ALAT
    use SystemData, only: nBasis!, iSpinSkip
    ! HACK - We use nBasisMax(2,3) here rather than iSpinSkip, as it appears
    !        to be more reliably set (see for example test H2O_RI)
    ! TODO: We need to sort this out so that they are consistent
    !       --> Talk to George/Alex to see what impact that might have?
    ! TODO: It would be nice to reduce the number of variants of sltcnd_...
    !       which are floating around.
    use HElem
    use UMatCache, only: GTID
    use IntegralsData, only: UMAT
    use OneEInts, only: GetTMatEl
    use Integrals, only: GetUMatEl
    use DetBitOps, only: count_open_orbs, FindBitExcitLevel
    use csf_data, only: csf_sort_det_block
    use timing
    implicit none

contains
    function sltcnd_compat (nI, nJ, IC) result (hel)

        ! Use the Slater-Condon Rules to evaluate the H-matrix element between
        ! two determinants, where the value of IC is already known.
        !
        ! In:  nI, nJ       - The determinants to evaluate
        !      IC           - The number of orbitals I,J differ by
        ! Ret: hel          - The H matrix element

        integer, intent(in) :: nI(nel), nJ(nel), IC
        type(HElement) :: hel

        integer :: ex(2,2)
        logical :: tParity

        select case (IC)
        case (0)
            ! The determinants are exactly the same
            hel = sltcnd_0 (nI)

        case (1)
            ! The determinants differ by only one orbital
            ex(1,1) = IC
            call GetExcitation (nI, nJ, nel, ex, tParity)
            hel = sltcnd_1 (nI, ex(:,1), tParity)

        case (2)
            ! The determinants differ by two orbitals
            ex(1,1) = IC
            call GetExcitation (nI, nJ, nel, ex, tParity)
            hel = sltcnd_2 (ex, tParity)

        case default
            ! The determinants differ by more than 2 orbitals
            hel = helement(0)
        end select
    end function sltcnd_compat


    type(HElement) function sltcnd_excit (nI, nJ, IC, ex, tParity)
        
        ! Use the Slater-Condon Rules to evaluate the H-matrix element between
        ! two determinants, where the excitation matrix is already known.
        !
        ! In:  nI, nJ       - The determinants to evaluate
        !      IC           - The number of orbitals I,J differ by
        !      ex           - The excitation matrix
        !      tParity      - The parity of the excitation
        ! Ret: sltcnd_excit - The H matrix element

        integer, intent(in) :: nI(nel), nJ(nel), IC
        integer, intent(in), optional :: ex(2,2)
        logical, intent(in), optional :: tParity
        character(*), parameter :: this_routine = 'sltcnd_excit'

        if (IC /= 0 .and. .not. present(tParity)) &
            call stop_all (this_routine, "ex and tParity must be provided to &
                          &sltcnd_excit for all IC /= 0")

        select case (IC)
        case (0)
            ! The determinants are exactly the same
            sltcnd_excit = sltcnd_0 (nI)

        case (1)
            ! The determnants differ by only one orbital
            sltcnd_excit = sltcnd_1 (nI, ex(:,1), tParity)

        case (2)
            ! The determinants differ by two orbitals
            sltcnd_excit = sltcnd_2 (ex, tParity)

        case default
            ! The determinants differ yb more than 2 orbitals
            sltcnd_excit%v = 0
        end select
    end function

    function sltcnd_knowIC (nI, nJ, iLutI, iLutJ, IC) result(hel)

        ! Use the Slater-Condon Rules to evaluate the H-matrix element between
        ! two determinants, where the value of IC and the bit representations
        ! are already known.
        !
        ! In:  nI, nJ       - The determinants to evaluate
        !      iLutI, iLutJ - Bit representations of I,J
        !      IC           - The number of orbitals I,J differ by
        ! Ret: hel          - The H matrix element

        integer, intent(in) :: nI(nel), nJ(nel)
        integer, intent(in) :: iLutI(0:NIfTot), iLutJ(0:NIfTot)
        integer, intent(in) :: IC
        type(HElement) :: hel
        integer :: ex(2,2)
        logical :: tSign

        select case (IC)
        case (0)
            ! The determinants are exactly the same
            hel = sltcnd_0 (nI)

        case (1)
            ! The determinants differ by only one orbital
            ex(1,1) = IC
            call GetBitExcitation (iLutI, iLutJ, ex, tSign)
            hel = sltcnd_1 (nI, Ex(:,1), tSign)

        case (2)
            ! The determinants differ by two orbitals
            ex(1,1) = IC
            call GetBitExcitation (iLutI, iLutJ, ex, tSign)
            hel = sltcnd_2 (ex, tSign)

        case default
            ! The determinants differ by more than two orbitals
            hel%v = 0
        endselect

    end function
    
    type(HElement) function sltcnd (nI, nJ, iLutI, iLutJ, ICret)
        
        ! Use the Slater-Condon Rules to evaluate the H matrix element between
        ! two determinants. Make no assumptions about ordering of orbitals.
        ! However, this is NOT to be passed CSFS - it is to evaluate the 
        ! component determinants.
        !
        ! In:  nI, nJ        - The determinants to evaluate
        !      iLutI, ilutJ  - Bit representations of above determinants
        ! Out: ICret         - Optionally return the IC value
        ! Ret: sltcnd        - The H matrix element

        integer, intent(in) :: nI(nel), nJ(nel)
        integer, intent(in) :: iLutI(0:NIfTot), iLutJ(0:NIfTot)
        integer, intent(out), optional :: ICret
        integer :: IC

        ! Get the excitation level
        IC = FindBitExcitLevel (iLutI, iLutJ)

        sltcnd = sltcnd_knowIC (nI, nJ, iLutI, iLutJ, IC)

        if (present(ICRet)) ICret = IC

    end function


    function sltcnd_0 (nI) result(hel)

        ! Calculate the HElement by the SlaterCondon Rules when the two
        ! determinants are the same (so we only need to specify one).

        integer, intent(in) :: nI(nel)
        type (HElement) :: hel, hel_sing, hel_doub, hel_tmp
        integer :: id(nel), ids, i, j, idN, idX

        ! Sum in the one electron integrals (KE --> TMAT)
        hel_sing = sum(GetTMATEl(nI, nI))

        ! Obtain the spatial rather than spin indices if required
        id = gtID(nI)

        ! Sum in the two electron contributions. Use max(id...) as we cannot
        ! guarantee that if j>i then nI(j)>nI(i).
        hel_doub = HElement(0)
        hel_tmp = HElement(0)
        do i=1,nel-1
            do j=i+1,nel
                idX = max(id(i), id(j))
                idN = min(id(i), id(j))
                hel_doub = hel_doub + GetUMATEl(nBasisMax, UMAT, ALAT, &
                                                nBasis, nBasisMax(2,3), G1, &
                                                idN, idX, idN, idX)
            enddo
        enddo
                
        ! Exchange contribution only considered if tExch set.
        ! This is only separated from the above loop to keep "if (tExch)" out
        ! of the tight loop for efficiency.
        if (tExch) then
            do i=1,nel-1
                do j=i+1,nel
                    ! Exchange contribution is zero if I,J are alpha/beta
                    if (G1(nI(i))%Ms == G1(nI(j))%Ms) then
                        idX = max(id(i), id(j))
                        idN = min(id(i), id(j))
                        hel_tmp = hel_tmp - GetUMATEl(nBasisMax, UMAT, ALAT, &
                                                      nBasis, nBasisMax(2,3),&
                                                      G1, idN, idX, idX, idN)
                    endif
                enddo
            enddo
        endif
        hel_doub = hel_doub + hel_tmp

        ! If we are scaling the coulomb interaction, do so here.
        hel = hel_sing + (hel_doub * HElement(FCOUL))
    end function sltcnd_0

    function sltcnd_1 (nI, ex, tSign) result(hel)

        ! Calculate the HElement by the Slater-Condon Rules when the two
        ! determinants differ by one orbital exactly.

        integer, intent(in) :: nI(nel), ex(2)
        logical, intent(in) :: tSign
        type (HElement) :: hel
        integer :: id_ex(2), id, i

        ! Obtain spatial rather than spin indices if required
        id_ex = gtID(ex)

        ! Sum in the diagonal terms (same in both dets)
        ! Coulomb term only included if Ms values of ex(1) and ex(2) are the
        ! same.
        hel = HElement(0)
        if (G1(ex(1))%Ms == G1(ex(2))%Ms) then
            do i=1,nel
                if (ex(1) /= nI(i)) then
                    id = gtID(nI(i))
                    hel = hel + GetUMATEl (nBasisMax, UMAT, ALAT, nBasis, &
                                           nBasisMax(2,3), G1, id_ex(1), id, &
                                           id_ex(2), id)
                endif
            enddo
        endif

        ! Exchange contribution is only considered if tExch set.
        ! This is only separated from the above loop to keep "if (tExch)" out
        ! of the tight loop for efficiency.
        if (tExch .and. G1(ex(1))%Ms == G1(ex(2))%Ms) then
            do i=1,nel
                if (ex(1) /= nI(i)) then
                    if (G1(ex(1))%Ms == G1(nI(i))%Ms) then
                        id = gtID(nI(i))
                        hel = hel - GetUMATEl (nBasisMax, UMAT, ALAT, nBasis,&
                                               nBasisMax(2,3), G1, id_ex(1),&
                                               id, id, id_ex(2))
                    endif
                endif
            enddo
        endif

        ! consider the non-diagonal part of the kinetic energy -
        ! <psi_a|T|psi_a'> where a, a' are the only basis fns that differ in
        ! nI, nJ
        hel = (hel*HElement(FCOUL)) + GetTMATEl(ex(1), ex(2))

        if (tSign) hel = -hel
    end function sltcnd_1
    
    function sltcnd_2 (ex, tSign) result (hel)

        ! Calculate the HElement by the Slater-Condon Rules when the two
        ! determinants differ by two orbitals exactly (the simplest case).

        integer, intent(in) :: ex(2,2)
        logical, intent(in) :: tSign
        type (HElement) :: hel
        integer :: id(2,2)

        ! Obtain spatial rather than spin indices if required
        id = gtID(ex)

        ! Only non-zero contributions if Ms preserved in each term (consider
        ! physical notation).
        if ( (G1(ex(1,1))%Ms == G1(ex(2,1))%Ms) .and. &
             (G1(ex(1,2))%Ms == G1(ex(2,2))%Ms) ) then
             hel = GetUMATEl (nBasisMax, UMAT, ALAT, nBasis, nBasisMax(2,3), &
                              G1, id(1,1), id(1,2), id(2,1), id(2,2))
        else
            hel = HElement(0)
        endif

        if ( (G1(ex(1,1))%Ms == G1(ex(2,2))%Ms) .and. &
             (G1(ex(1,2))%Ms == G1(Ex(2,1))%Ms) ) then
             hel = hel - GetUMATEl (nBasismax, UMAT, ALAT, nBasis, &
                                    nBasisMax(2,3), G1, id(1,1), id(1,2), &
                                    id(2,2), id(2,1))
        endif

        if (tSign) hel = -hel
    end function sltcnd_2
end module
