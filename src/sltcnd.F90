#include "macros.h"

module sltcnd_mod
    use SystemData, only: nel, nBasisMax, tExch, G1, ALAT, tReltvy
    use SystemData, only: nBasis!, iSpinSkip
    ! HACK - We use nBasisMax(2,3) here rather than iSpinSkip, as it appears
    !        to be more reliably set (see for example test H2O_RI)
    ! TODO: We need to sort this out so that they are consistent
    !       --> Talk to George/Alex to see what impact that might have?
    ! TODO: It would be nice to reduce the number of variants of sltcnd_...
    !       which are floating around.
    use constants, only: dp,n_int
    use UMatCache, only: GTID
    use IntegralsData, only: UMAT
    use OneEInts, only: GetTMatEl
    use Integrals_neci, only: get_umat_el
    use DetBitOps, only: count_open_orbs, FindBitExcitLevel
    use csf_data, only: csf_sort_det_block
    use timing_neci
    use bit_reps, only: NIfTot
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
        HElement_t(dp) :: hel

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
            hel = (0)
        end select
    end function sltcnd_compat


    HElement_t(dp) function sltcnd_excit (nI, IC, ex, tParity)

        ! Use the Slater-Condon Rules to evaluate the H-matrix element between
        ! two determinants, where the excitation matrix is already known.
        !
        ! In:  nI, nJ       - The determinants to evaluate
        !      IC           - The number of orbitals I,J differ by
        !      ex           - The excitation matrix
        !      tParity      - The parity of the excitation
        ! Ret: sltcnd_excit - The H matrix element

        integer, intent(in) :: nI(nel), IC
        integer, intent(in), optional :: ex(2,2)
        logical, intent(in), optional :: tParity
        character(*), parameter :: this_routine = 'sltcnd_excit'

        if (IC /= 0 .and. .not. (present(ex) .and. present(tParity))) &
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
            sltcnd_excit = 0
        end select
    end function

    function sltcnd_knowIC (nI, iLutI, iLutJ, IC) result(hel)

        ! Use the Slater-Condon Rules to evaluate the H-matrix element between
        ! two determinants, where the value of IC and the bit representations
        ! are already known.
        !
        ! In:  nI, nJ       - The determinants to evaluate
        !      iLutI, iLutJ - Bit representations of I,J
        !      IC           - The number of orbitals I,J differ by
        ! Ret: hel          - The H matrix element

        integer, intent(in) :: nI(nel)
        integer(kind=n_int), intent(in) :: iLutI(0:NIfTot), iLutJ(0:NIfTot)
        integer, intent(in) :: IC
        HElement_t(dp) :: hel
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
            hel = 0
        endselect

    end function

    HElement_t(dp) function sltcnd (nI, iLutI, iLutJ, ICret)

        ! Use the Slater-Condon Rules to evaluate the H matrix element between
        ! two determinants. Make no assumptions about ordering of orbitals.
        ! However, this is NOT to be passed CSFS - it is to evaluate the
        ! component determinants.
        !
        ! In:  nI, nJ        - The determinants to evaluate
        !      iLutI, ilutJ  - Bit representations of above determinants
        ! Out: ICret         - Optionally return the IC value
        ! Ret: sltcnd        - The H matrix element

        integer, intent(in) :: nI(nel)
        integer(kind=n_int), intent(in) :: iLutI(0:NIfTot), iLutJ(0:NIfTot)
        integer, intent(out), optional :: ICret
        integer :: IC

        ! Get the excitation level
        IC = FindBitExcitLevel (iLutI, iLutJ)

        sltcnd = sltcnd_knowIC (nI, iLutI, iLutJ, IC)

        if (present(ICRet)) ICret = IC

    end function

    function CalcFockOrbEnergy (Orb,HFDet) result(hel)
        ! This calculates the orbital fock energy from
        ! the one- and two- electron integrals. This
        ! requires a knowledge of the HF determinant.
        !In: Orbital (Spin orbital notation)
        !In: HFDet (HF Determinant)
        integer, intent(in) :: HFDet(nel),Orb
        integer :: idHF(NEl),idOrb,j,idN
        HElement_t(dp) :: hel_sing,hel

        !GetTMATEl works with spin orbitals
        hel_sing = GetTMATEl(Orb,Orb)

        ! Obtain the spatial rather than spin indices if required
        idOrb = gtID(Orb)
        idHF = gtID(HFDet)

        ! Sum in the two electron contributions.
        hel = (0)
        do j=1,nel
            idN = idHF(j)
            hel = hel + get_umat_el (idOrb, idN, idOrb, idN)
        enddo

        ! Exchange contribution only considered if tExch set.
        ! This is only separated from the above loop to keep "if (tExch)" out
        ! of the tight loop for efficiency.
        do j=1,nel
            ! Exchange contribution is zero if I,J are alpha/beta
            if (tReltvy.or.(G1(Orb)%Ms == G1(HFDet(j))%Ms)) then
                idN = idHF(j)
                hel = hel - get_umat_el (idOrb, idN, idN, idOrb)
            endif
        enddo
        hel = hel + hel_sing

    end function CalcFockOrbEnergy

    function SumFock (nI,HFDet) result(hel)

        ! This just calculates the sum of the Fock energies
        ! by considering the one-electron integrals and
        ! the double-counting contribution
        ! to the diagonal matrix elements. This is subtracted from
        ! the sum of the fock energies to calculate diagonal
        ! matrix elements, or added to the sum of the 1-electron
        ! integrals. The HF determinant needs to be supplied.

        integer , intent(in) :: nI(nel),HFDet(nel)
        HElement_t(dp) :: hel,hel_doub,hel_tmp,hel_sing
        integer :: i,j,idN,idX,id(nel),idHF(NEl)

        !Obtain the 1e terms
        hel_sing = sum(GetTMATEl(nI, nI))

        ! Obtain the spatial rather than spin indices if required
        id = gtID(nI)
        idHF = gtID(HFDet)

        ! Sum in the two electron contributions. Use max(id...) as we cannot
        ! guarantee that if j>i then nI(j)>nI(i).
        hel_doub = (0)
        hel_tmp = (0)
        do i=1,nel
            do j=1,nel
                idX = id(i)
                idN = idHF(j)
                hel_doub = hel_doub + get_umat_el (idX, idN, idX, idN)
            enddo
        enddo

        ! Exchange contribution only considered if tExch set.
        ! This is only separated from the above loop to keep "if (tExch)" out
        ! of the tight loop for efficiency.
        if (tExch) then
            do i=1,nel
                do j=1,nel
                    ! Exchange contribution is zero if I,J are alpha/beta
                    if (tReltvy.or.(G1(nI(i))%Ms == G1(HFDet(j))%Ms)) then
                        idX = id(i)
                        idN = idHF(j)
                        hel_tmp = hel_tmp - get_umat_el (idX, idN, idN, idX)
                    endif
                enddo
            enddo
        endif
        hel = hel_doub + hel_tmp + hel_sing

    end function SumFock

    function sltcnd_0 (nI) result(hel)

        ! Calculate the  by the SlaterCondon Rules when the two
        ! determinants are the same (so we only need to specify one).

        integer, intent(in) :: nI(nel)
        HElement_t(dp) :: hel, hel_sing, hel_doub, hel_tmp
        integer :: id(nel), i, j, idN, idX

        ! Sum in the one electron integrals (KE --> TMAT)
        hel_sing = sum(GetTMATEl(nI, nI))

        ! Obtain the spatial rather than spin indices if required
        id = gtID(nI)
        !write(6,*) "****",id(:)

        ! Sum in the two electron contributions. Use max(id...) as we cannot
        ! guarantee that if j>i then nI(j)>nI(i).
        hel_doub = (0)
        hel_tmp = (0)
        do i=1,nel-1
            do j=i+1,nel
                idX = max(id(i), id(j))
                idN = min(id(i), id(j))
                hel_doub = hel_doub + get_umat_el (idN, idX, idN, idX)
                !write(6,*) idN,idX,idN,idX,get_umat_el (idN,idX,idN,idX)
            enddo
        enddo

        ! Exchange contribution only considered if tExch set.
        ! This is only separated from the above loop to keep "if (tExch)" out
        ! of the tight loop for efficiency.
        if (tExch) then
            do i=1,nel-1
                do j=i+1,nel
                    ! Exchange contribution is zero if I,J are alpha/beta
                    if ((G1(nI(i))%Ms == G1(nI(j))%Ms).or.tReltvy) then
                        idX = max(id(i), id(j))
                        idN = min(id(i), id(j))
                        hel_tmp = hel_tmp - get_umat_el (idN, idX, idX, idN)
                        !write(6,*) idN,idX,idX,idN,get_umat_el (idN,idX,idX,idN)
                    endif
                enddo
            enddo
        endif
        hel = hel_doub + hel_tmp + hel_sing

    end function sltcnd_0

    function sltcnd_1 (nI, ex, tSign) result(hel)

        ! Calculate the  by the Slater-Condon Rules when the two
        ! determinants differ by one orbital exactly.

        integer, intent(in) :: nI(nel), ex(2)
        logical, intent(in) :: tSign
        HElement_t(dp) :: hel
        integer :: id_ex(2), id, i

        ! Obtain spatial rather than spin indices if required
        id_ex = gtID(ex)

        ! Sum in the diagonal terms (same in both dets)
        ! Coulomb term only included if Ms values of ex(1) and ex(2) are the
        ! same.
        hel = (0)
        if (tReltvy.or.(G1(ex(1))%Ms == G1(ex(2))%Ms)) then
            do i=1,nel
                if (ex(1) /= nI(i)) then
                    id = gtID(nI(i))
                    hel = hel + get_umat_el (id_ex(1), id, id_ex(2), id)
                endif
            enddo
        endif

        ! Exchange contribution is only considered if tExch set.
        ! This is only separated from the above loop to keep "if (tExch)" out
        ! of the tight loop for efficiency.
        if (tExch .and. ((G1(ex(1))%Ms == G1(ex(2))%Ms).or.tReltvy)) then
            do i=1,nel
                if (ex(1) /= nI(i)) then
                    if (tReltvy.or.(G1(ex(1))%Ms == G1(nI(i))%Ms)) then
                        id = gtID(nI(i))
                        hel = hel - get_umat_el (id_ex(1), id, id, id_ex(2))
                    endif
                endif
            enddo
        endif

        ! consider the non-diagonal part of the kinetic energy -
        ! <psi_a|T|psi_a'> where a, a' are the only basis fns that differ in
        ! nI, nJ
        hel = hel + GetTMATEl(ex(1), ex(2))

        if (tSign) hel = -hel
    end function sltcnd_1

    function sltcnd_2 (ex, tSign) result (hel)

        ! Calculate the  by the Slater-Condon Rules when the two
        ! determinants differ by two orbitals exactly (the simplest case).

        integer, intent(in) :: ex(2,2)
        logical, intent(in) :: tSign
        HElement_t(dp) :: hel
        integer :: id(2,2)

        ! Obtain spatial rather than spin indices if required
        id = gtID(ex)

        ! Only non-zero contributions if Ms preserved in each term (consider
        ! physical notation).
!>>>!        write(6, '("---> ")', advance='no')
        if ( tReltvy.or.((G1(ex(1,1))%Ms == G1(ex(2,1))%Ms) .and. &
             (G1(ex(1,2))%Ms == G1(ex(2,2))%Ms)) ) then
             hel = get_umat_el (id(1,1), id(1,2), id(2,1), id(2,2))
        else
            hel = (0)
        endif
!>>>!        write(6,'(f10.6)', advance='no') hel

        if ( tReltvy.or.((G1(ex(1,1))%Ms == G1(ex(2,2))%Ms) .and. &
             (G1(ex(1,2))%Ms == G1(Ex(2,1))%Ms)) ) then
             hel = hel - get_umat_el (id(1,1), id(1,2), id(2,2), id(2,1))
!>>>!            write(6,'(f10.6)', advance='no') - get_umat_el (id(1,1), id(1,2), id(2,2), id(2,1))
!>>>!        else
!>>>!            write(6,'(f10.6)', advance='no') 0.0
        endif
!>>>!        write(6,*) hel

        if (tSign) hel = -hel
    end function sltcnd_2
end module
