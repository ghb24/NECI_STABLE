module sltcnd_csf_mod
    use SystemData, only: nel, nBasisMax, tExch, FCOUL, NIfTot, G1, ALAT
    use SystemData, only: nBasis, iSpinSkip
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
    type(HElement) function sltcnd_csf (nI, nJ, iLutI, iLutJ)
        
        ! Use the Slater-Condon Rules to evaluate the H matrix element between
        ! two determinants. Assume CSF ordering of orbitals (closed pairs
        ! followed by open shell electrons). However, this is NOT to be passed
        ! CSFS - it is to evaluate the component determinants.
        !
        ! In:  nI, nJ        - The determinants to evaluate
        !      iLutI, ilutJ  - Bit representations of above determinants
        ! Ret: sltcnd_csf    - The H matrix element

        integer, intent(in) :: nI(nel), nJ(nel)
        integer, intent(in) :: iLutI(0:NIfTot), iLutJ(0:NIfTot)
        integer :: IC, ex(2,2)
        logical :: tSign
        
        ! Get the excitation level
        IC = FindBitExcitLevel (iLutI, iLutJ)

        select case (IC)
        case (0)
            ! The determinants are exactly the same
            sltcnd_csf = sltcnd_csf_0 (nI)

        case (1)
            ! The determinants differ by only one orbital
            ! TODO: For speed, can do a no-tSign version...
            ex(1,1) = IC
            call GetBitExcitation (iLutI, iLutJ, ex, tSign)
            sltcnd_csf = sltcnd_csf_1 (nI, Ex(:,1), tSign)

        case (2)
            ! The determinants differ by two orbitals
            ex(1,1) = IC
            call GetBitExcitation (iLutI, iLutJ, ex, tSign)
            sltcnd_csf = sltcnd_csf_2 (ex, tSign)

        case default
            ! The determinants differ by more than two orbitals
            sltcnd_csf%v = 0
        endselect

    end function


    function sltcnd_csf_0 (nI) result(hel)

        ! Calculate the HElement by the SlaterCondon Rules when the two
        ! determinants are the same (so we only need to specify one).

        integer, intent(in) :: nI(nel)
        type (HElement) :: hel, hel_sing, hel_doub, hel_tmp
        integer :: id(nel), ids, i, j, idN, idX

        ! Sum in the one electron integrals (KE --> TMAT)
        hel_sing = HElement(0)
        do i=1,nel
            hel_sing = hel_sing + GetTMATEl(nI(i), nI(i))
        enddo

        ! Obtain the spatial rather than spin indices if required
        id = gtID(nI)

        ! Sum in the two electron contributions. Use max(id...) as we cannot
        ! guarantee that if j>i then nI(j)>nI(i).
        hel_doub = HElement(0)
        hel_tmp = HElement(0)
        do i=1,nel-1
            do j=i,nel
                idX = max(id(i), id(j))
                idN = min(id(i), id(j))
                hel_doub = hel_doub + GetUMATEl(nBasisMax, UMAT, ALAT, &
                                                nBasis, iSpinSkip, G1, idN, &
                                                idX, idN, idX)
                
                ! If are not considering the exchange contribution, or if I,J
                ! are alpha/beta (ie exchange == 0) then don't continue
                ! TODO: Remove ability to turn off exchange (tExch)?
                ids = G1(nI(i))%Ms * G1(nI(j))%Ms
                if (tExch .and. ids > 0) then
                    hel_tmp = hel_tmp - GetUMATEl(nBasisMax, UMAT, ALAT, &
                                                  nBasis, iSpinSkip, G1, idN,&
                                                  idX, idX, idN)
                endif
            enddo
        enddo
        hel_doub = hel_doub + hel_tmp

        ! If we are scaling the coulomb interaction, do so here.
        ! TODO: Should we remove ability to use FCOUL
        hel = hel_sing + (hel_doub * HElement(FCOUL))
    end function sltcnd_csf_0

    function sltcnd_csf_1 (nI, ex, tSign) result(hel)

        ! Calculate the HElement by the Slater-Condon Rules when the two
        ! determinants differ by one orbital exactly.

        integer, intent(in) :: nI(nel), ex(2)
        logical, intent(in) :: tSign
        type (HElement) :: hel
        integer :: id_ex(2), id, i
        logical :: bEqualMs

        ! Obtain spatial rather than spin indices if required
        id_ex = gtID(ex)

        ! Sum in the diagonal terms (same in both dets)
        hel = HElement(0)
        bEqualMs = G1(ex(1))%Ms == G1(ex(2))%Ms
        do i=1,nel
            if (ex(1) /= nI(i)) then
                id = gtID(nI(i))
                if (bEqualMs) then
                    hel = hel + GetUMATEl (nBasisMax, UMAT, ALAT, nBasis, &
                                           iSpinSkip, G1, id_ex(1), id, &
                                           id_ex(2), id)
                endif

                ! TODO: Remove tExch
                if (tExch .and. (G1(ex(1))%Ms == G1(nI(i))%Ms) .and. &
                                (G1(ex(2))%Ms == G1(nI(i))%Ms) ) then
                    hel = hel - GetUMATEl (nBasisMax, UMAT, ALAT, nBasis, &
                                           iSpinSkip, G1, id_ex(1), id, &
                                           id, id_ex(2))
                endif
            endif
        enddo

        ! consider the non-diagonal part of the kinetic energy -
        ! <psi_a|T|psi_a'> where a, a' are the only basis fns that differ in
        ! nI, nJ
        ! TODO: Remove FCOUL?
        hel = (hel*HElement(FCOUL)) + GetTMATEl(ex(1), ex(2))

        if (tSign) hel = -hel
    end function sltcnd_csf_1
    
    function sltcnd_csf_2 (ex, tSign) result (hel)

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
             hel = GetUMATEl (nBasisMax, UMAT, ALAT, nBasis, iSpinSkip, G1, &
                              id(1,1), id(1,2), id(2,1), id(2,2))
        else
            hel = HElement(0)
        endif

        if ( (G1(ex(1,1))%Ms == G1(ex(2,2))%Ms) .and. &
             (G1(ex(1,2))%Ms == G1(Ex(2,1))%Ms) ) then
             hel = hel - GetUMATEl (nBasismax, UMAT, ALAT, nBasis, iSpinSkip,&
                                    G1, id(1,1), id(1,2), id(2,2), id(2,1))
        endif

        if (tSign) hel = -hel
    end function sltcnd_csf_2
end module
