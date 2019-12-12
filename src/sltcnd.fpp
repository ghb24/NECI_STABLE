#include "macros.h"
#:include "macros.fpph"
#:set ExcitationTypes = ['NoExc_t', 'SingleExc_t', 'DoubleExc_t', 'TripleExc_t']

module sltcnd_mod
    use SystemData, only: nel, nBasisMax, tExch, G1, ALAT, tReltvy, t_3_body_excits, &
                          t_mol_3_body, t_ueg_3_body, tContact
    use SystemData, only: nBasis!, iSpinSkip
    ! HACK - We use nBasisMax(2,3) here rather than iSpinSkip, as it appears
    !        to be more reliably set (see for example test H2O_RI)
    ! TODO: We need to sort this out so that they are consistent
    !       --> Talk to George/Alex to see what impact that might have?
    use constants, only: dp, n_int, maxExcit
    use UMatCache, only: GTID, UMatInd
    use IntegralsData, only: UMAT
    use OneEInts, only: GetTMatEl, TMat2D
    use procedure_pointers, only: get_umat_el
    use excitation_types, only: NoExc_t, SingleExc_t, DoubleExc_t, TripleExc_t, UNKNOWN
    use DetBitOps, only: count_open_orbs, FindBitExcitLevel
    use csf_data, only: csf_sort_det_block
    use timing_neci
    use bit_reps, only: NIfTot
    use LMat_mod, only: get_lmat_el, get_lmat_el_ua
    use gen_coul_ueg_mod, only: get_contact_umat_el_3b_sp, get_contact_umat_el_3b_sap
    use tc_three_body_data, only: tDampKMat, tSpinCorrelator
    use kMatProjE, only: kMatParSpinCorrection, kMatOppSpinCorrection, spinKMatContrib, &
                         kMatAA
    implicit none
    private
    public :: initSltCndPtr, sltcnd_compat, sltcnd, CalcFockOrbEnergy, &
              sltcnd_2_kernel, &
              sumfock, sltcnd_knowIC, &
              sltcnd_excit, sltcnd_excit_old

!> Evaluate the H-matrix element using the Slater-Condon rules.
!> Generic function that accepts arguments of ExcitationTypes,
!> or the old-style integers.
    interface sltcnd_excit
    #:for excitation_t in ExcitationTypes
        module procedure sltcnd_excit_${excitation_t}$
    #:endfor
    end interface

    abstract interface
        function sltcnd_0_t(nI, exc) result(hel)
            import :: dp, nel, NoExc_t
            integer, intent(in) :: nI(nel)
            type(NoExc_t), intent(in) :: exc
            HElement_t(dp) :: hel
        end function sltcnd_0_t

        function sltcnd_1_t(nI, exc, tSign) result(hel)
            import :: dp, nel, SingleExc_t
            integer, intent(in) :: nI(nel)
            type(SingleExc_t), intent(in) :: exc
            logical, intent(in) :: tSign
            HElement_t(dp) :: hel
        end function sltcnd_1_t

        function sltcnd_2_t(nI, ex, tSign) result(hel)
            import :: dp, nel, DoubleExc_t
            integer, intent(in) :: nI(nel)
            type(DoubleExc_t), intent(in) :: ex
            logical, intent(in) :: tSign
            HElement_t(dp) :: hel
        end function sltcnd_2_t

        function sltcnd_3_t(ex,tSign) result(hel)
            import :: dp, nel
            integer, intent(in) :: ex(2,3)
            logical, intent(in) :: tSign
            HElement_t(dp) :: hel
        end function sltcnd_3_t
    end interface

    procedure(sltcnd_0_t), pointer :: sltcnd_0
    procedure(sltcnd_1_t), pointer :: sltcnd_1
    procedure(sltcnd_2_t), pointer :: sltcnd_2
    procedure(sltcnd_3_t), pointer :: sltcnd_3

contains

    subroutine initSltCndPtr()
        use SystemData, only: tSmallBasisForThreeBody
        implicit none

        if (TContact) then

            if (t_mol_3_body &
                .or. t_ueg_3_body .and. nel > 2 .and. tSmallBasisForThreeBody) then

                sltcnd_0 => sltcnd_0_tc_ua
                sltcnd_1 => sltcnd_1_tc_ua
                sltcnd_2 => sltcnd_2_tc_ua
                sltcnd_3 => sltcnd_3_tc_ua
            else
                sltcnd_0 => sltcnd_0_base_ua
                sltcnd_1 => sltcnd_1_base_ua
                sltcnd_2 => sltcnd_2_base_ua
                sltcnd_3 => sltcnd_3_base
            end if
        else
            ! six-index integrals are only used for three and more
            ! electrons
            if (t_mol_3_body .or. t_ueg_3_body .and. nel > 2) then
                sltcnd_0 => sltcnd_0_tc
                sltcnd_1 => sltcnd_1_tc
                sltcnd_2 => sltcnd_2_tc
                sltcnd_3 => sltcnd_3_tc
            else
                sltcnd_0 => sltcnd_0_base
                sltcnd_1 => sltcnd_1_base
                sltcnd_2 => sltcnd_2_base
                sltcnd_3 => sltcnd_3_base
            end if

        endif
    end subroutine initSltCndPtr

    function sltcnd_compat(nI, nJ, IC) result(hel)

        ! Use the Slater-Condon Rules to evaluate the H-matrix element between
        ! two determinants, where the value of IC is already known.
        !
        ! In:  nI, nJ       - The determinants to evaluate
        !      IC           - The number of orbitals I,J differ by
        ! Ret: hel          - The H matrix element

        integer, intent(in) :: nI(nel), nJ(nel), IC
        HElement_t(dp) :: hel
#ifdef DEBUG_
        character(*), parameter :: this_routine = "sltcnd_compat"
#endif

        integer :: ex(2, maxExcit)
        logical :: tParity

        select case (IC)
        case (0)
            ! The determinants are exactly the same
            hel = sltcnd_excit(nI, NoExc_t())

        case (1)
            ! The determinants differ by only one orbital
            ex(1, 1) = IC
            call GetExcitation(nI, nJ, nel, ex, tParity)
            hel = sltcnd_excit(nI, SingleExc_t(ex(:, 1)), tParity)
        case (2)
            ! The determinants differ by two orbitals
            ex(1, 1) = IC
            call GetExcitation(nI, nJ, nel, ex, tParity)
            hel = sltcnd_excit(nI, DoubleExc_t(ex(:, :2)), tParity)
        case (3)
            ! The determinants differ by three orbitals
            ex(1, 1) = IC
            call GetExcitation(nI, nJ, nel, ex, tParity)
            hel = sltcnd_excit(nI, TripleExc_t(ex(:, :3)), tParity)

        case default
            ! The determinants differ by more than 3 orbitals
            ASSERT(.not. t_3_body_excits)
            hel = h_cast(0.0_dp)
        end select
    end function sltcnd_compat

    function sltcnd_excit_old(nI, IC, ex, tParity) result(hel)

        ! Use the Slater-Condon Rules to evaluate the H-matrix element between
        ! two determinants, where the excitation matrix is already known.
        !
        ! In:  nI, nJ       - The determinants to evaluate
        !      IC           - The number of orbitals I,J differ by
        !      ex           - The excitation matrix
        !      tParity      - The parity of the excitation
        ! Ret: sltcnd_excit - The H matrix element

        integer, intent(in) :: nI(nel), IC
        integer, intent(in), optional :: ex(2, ic)
        logical, intent(in), optional :: tParity
        HElement_t(dp) :: hel
        character(*), parameter :: this_routine = 'sltcnd_excit_old'

        if (IC /= 0 .and. .not. (present(ex) .and. present(tParity))) &
            call stop_all(this_routine, "ex and tParity must be provided to &
                          &sltcnd_excit for all IC /= 0")
        select case (IC)
        case (0)
            ! The determinants are exactly the same
            hel = sltcnd_excit(nI, NoExc_t())
        case (1)
            hel = sltcnd_excit(nI, SingleExc_t(ex(:, 1)), tParity)
        case (2)
            ! The determinants differ by two orbitals
            hel = sltcnd_excit(nI, DoubleExc_t(ex(:, :2)), tParity)
        case (3)
            hel = sltcnd_excit(nI, TripleExc_t(ex(:, :3)), tParity)
        case default
            ! The determinants differ by more than 2 orbitals
            ASSERT(.not. t_3_body_excits)
            hel = h_cast(0.0_dp)
        end select
    end function

    function sltcnd_knowIC(nI, iLutI, iLutJ, IC) result(hel)

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
        integer :: ex(2, maxExcit)
        logical :: tSign
#ifdef DEBUG_
        character(*), parameter :: this_routine = "sltcnd_knowIC"
#endif

        select case (IC)
        case (0)
            ! The determinants are exactly the same
            hel = sltcnd_excit(nI, NoExc_t())

        case (1)
            ! The determinants differ by only one orbital
            ex(1, 1) = IC
            call GetBitExcitation(iLutI, iLutJ, ex, tSign)
            hel = sltcnd_excit(nI, SingleExc_t(ex(:, 1)), tSign)

        case (2)
            ! The determinants differ by two orbitals
            ex(1, 1) = IC
            call GetBitExcitation(iLutI, iLutJ, ex, tSign)
            hel = sltcnd_excit(nI, DoubleExc_t(ex(:, :2)), tSign)
        case (3)
            ex(1, 1) = IC
            call GetBitExcitation(iLutI, iLutJ, ex, tSign)
            hel = sltcnd_excit(nI, TripleExc_t(ex(:, :3)), tSign)

        case default
            ! The determinants differ by more than two orbitals
            hel = h_cast(0.0_dp)
        endselect

    end function

    HElement_t(dp) function sltcnd(nI, iLutI, iLutJ, ICret)

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
        IC = FindBitExcitLevel(iLutI, iLutJ)

        sltcnd = sltcnd_knowIC(nI, iLutI, iLutJ, IC)

        if (present(ICRet)) ICret = IC

    end function

    function CalcFockOrbEnergy(Orb, HFDet) result(hel)
        ! This calculates the orbital fock energy from
        ! the one- and two- electron integrals. This
        ! requires a knowledge of the HF determinant.
        !In: Orbital (Spin orbital notation)
        !In: HFDet (HF Determinant)
        integer, intent(in) :: HFDet(nel), Orb
        integer :: idHF(NEl), idOrb, j, idN
        HElement_t(dp) :: hel_sing, hel

        !GetTMATEl works with spin orbitals
        hel_sing = GetTMATEl(Orb, Orb)

        ! Obtain the spatial rather than spin indices if required
        idOrb = gtID(Orb)
        idHF = gtID(HFDet)

        ! Sum in the two electron contributions.
        hel = (0)
        do j = 1, nel
            idN = idHF(j)
            hel = hel + get_umat_el(idOrb, idN, idOrb, idN)
        enddo

        ! Exchange contribution only considered if tExch set.
        ! This is only separated from the above loop to keep "if (tExch)" out
        ! of the tight loop for efficiency.
        do j = 1, nel
            ! Exchange contribution is zero if I,J are alpha/beta
            if (tReltvy .or. (G1(Orb)%Ms == G1(HFDet(j))%Ms)) then
                idN = idHF(j)
                hel = hel - get_umat_el(idOrb, idN, idN, idOrb)
            endif
        enddo
        hel = hel + hel_sing

    end function CalcFockOrbEnergy

    function SumFock(nI, HFDet) result(hel)

        ! This just calculates the sum of the Fock energies
        ! by considering the one-electron integrals and
        ! the double-counting contribution
        ! to the diagonal matrix elements. This is subtracted from
        ! the sum of the fock energies to calculate diagonal
        ! matrix elements, or added to the sum of the 1-electron
        ! integrals. The HF determinant needs to be supplied.

        integer, intent(in) :: nI(nel), HFDet(nel)
        HElement_t(dp) :: hel, hel_doub, hel_tmp, hel_sing
        integer :: i, j, idN, idX, id(nel), idHF(NEl)

        !Obtain the 1e terms
        hel_sing = sum(GetTMATEl(nI, nI))

        ! Obtain the spatial rather than spin indices if required
        id = gtID(nI)
        idHF = gtID(HFDet)

        ! Sum in the two electron contributions. Use max(id...) as we cannot
        ! guarantee that if j>i then nI(j)>nI(i).
        hel_doub = (0)
        hel_tmp = (0)
        do i = 1, nel
            do j = 1, nel
                idX = id(i)
                idN = idHF(j)
                hel_doub = hel_doub + get_umat_el(idX, idN, idX, idN)
            enddo
        enddo

        ! Exchange contribution only considered if tExch set.
        ! This is only separated from the above loop to keep "if (tExch)" out
        ! of the tight loop for efficiency.
        if (tExch) then
            do i = 1, nel
                do j = 1, nel
                    ! Exchange contribution is zero if I,J are alpha/beta
                    if (tReltvy .or. (G1(nI(i))%Ms == G1(HFDet(j))%Ms)) then
                        idX = id(i)
                        idN = idHF(j)
                        hel_tmp = hel_tmp - get_umat_el(idX, idN, idN, idX)
                    endif
                enddo
            enddo
        endif
        hel = hel_doub + hel_tmp + hel_sing

    end function SumFock

    function sltcnd_0_base(nI, exc) result(hel)
        ! Calculate the  by the SlaterCondon Rules when the two
        ! determinants are the same (so we only need to specify one).
        integer, intent(in) :: nI(nel)
        type(NoExc_t), intent(in) :: exc
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
        do i = 1, nel - 1
            do j = i + 1, nel
                hel_doub = hel_doub + get_umat_el(id(i), id(j), id(i), id(j))
                if (tSpinCorrelator) then
                    hel_doub = hel_doub &
       + spinKMatContrib(id(i), id(j), id(i), id(j), G1(nI(i))%MS, G1(nI(j))%MS)
                endif
                if (tDampKMat) then
                    ! same-spin correction: a factor of 0.5
                    if (G1(nI(i))%Ms == G1(nI(j))%Ms) then
         hel_doub = hel_doub + kMatParSpinCorrection(id(i), id(j), id(i), id(j))
                        ! opposite-spin correction (if different orbs)
                    else if (id(i) /= id(j)) then
         hel_doub = hel_doub + kMatOppSpinCorrection(id(i), id(j), id(i), id(j))
                    endif
                endif
            enddo
        enddo

        ! Exchange contribution only considered if tExch set.
        ! This is only separated from the above loop to keep "if (tExch)" out
        ! of the tight loop for efficiency.
        if (tExch) then
            do i = 1, nel - 1
                do j = i + 1, nel
                    ! Exchange contribution is zero if I,J are alpha/beta
                    if ((G1(nI(i))%Ms == G1(nI(j))%Ms) .or. tReltvy) then
                     hel_tmp = hel_tmp - get_umat_el(id(i), id(j), id(j), id(i))
                        if (tSpinCorrelator) &
                            ! here, i and j always have the same spin
              hel_tmp = hel_tmp + kMatAA%exchElement(id(i), id(j), id(i), id(j))
                        if (tDampKMat) &
           hel_tmp = hel_tmp - kMatParSpinCorrection(id(i), id(j), id(j), id(i))
                    endif
                enddo
            enddo
        endif
        hel = hel_doub + hel_tmp + hel_sing

    end function sltcnd_0_base

    function sltcnd_1_base(nI, ex, tSign) result(hel)
        ! Calculate the  by the Slater-Condon Rules when the two
        ! determinants differ by one orbital exactly.
        integer, intent(in) :: nI(nel)
        type(SingleExc_t), intent(in) :: ex
        logical, intent(in) :: tSign
        HElement_t(dp) :: hel

        ! Sum in the diagonal terms (same in both dets)
        ! Coulomb term only included if Ms values of ex(1) and ex(2) are the
        ! same.

        hel = sltcnd_1_kernel(nI, ex)
        if (tSign) hel = -hel
    end function sltcnd_1_base

    function sltcnd_1_kernel(nI, ex) result(hel)
        integer, intent(in) :: nI(nel)
        type(SingleExc_t), intent(in) :: ex
        HElement_t(dp) :: hel
        integer :: i, id, id_ex(2)

        ! Obtain spatial rather than spin indices if required
        id_ex = gtID(ex%val)

        hel = (0)
        if (tReltvy .or. (G1(ex%val(1))%Ms == G1(ex%val(2))%Ms)) then
            do i = 1, nel
                if (ex%val(1) /= nI(i)) then
                    id = gtID(nI(i))
                    hel = hel + get_umat_el(id_ex(1), id, id_ex(2), id)
                    if (tSpinCorrelator) then
                        hel = hel &
       + spinKMatContrib(id_ex(1), id, id_ex(2), id, G1(ex%val(1))%MS, G1(nI(i))%MS)
                    endif
                    if (tDampKMat) then
                        if (G1(ex%val(1))%Ms == G1(nI(i))%Ms) then
                   hel = hel + kMatParSpinCorrection(id_ex(1), id, id_ex(2), id)
                        else
                            ! opposite spin correction only for different orbs, which is given here
                   hel = hel + kMatOppSpinCorrection(id_ex(1), id, id_ex(2), id)
                        endif
                    endif
                endif
            enddo
        endif
        ! Exchange contribution is only considered if tExch set.
        ! This is only separated from the above loop to keep "if (tExch)" out
        ! of the tight loop for efficiency.
        if (tExch .and. ((G1(ex%val(1))%Ms == G1(ex%val(2))%Ms) .or. tReltvy)) then
            do i = 1, nel
                if (ex%val(1) /= nI(i)) then
                    if (tReltvy .or. (G1(ex%val(1))%Ms == G1(nI(i))%Ms)) then
                        id = gtID(nI(i))
                        hel = hel - get_umat_el(id_ex(1), id, id, id_ex(2))
                        if (tSpinCorrelator) &
                      hel = hel + kMatAA%exchElement(id_ex(1), id, id_ex(2), id)
                        if (tDampKMat) &
                            ! only parallel spin excits can have exchange terms entering in the one-body
                            ! matrix element -> no need for opposite spin correction
                   hel = hel - kMatParSpinCorrection(id_ex(1), id, id, id_ex(2))
                    endif
                endif
            enddo
        endif
        ! consider the non-diagonal part of the kinetic energy -
        ! <psi_a|T|psi_a'> where a, a' are the only basis fns that differ in
        ! nI, nJ
        hel = hel + GetTMATEl(ex%val(1), ex%val(2))
    end function sltcnd_1_kernel

    function sltcnd_2_base(nI, exc, tSign) result(hel)
        ! Calculate the  by the Slater-Condon Rules when the two
        ! determinants differ by two orbitals exactly (the simplest case).
        integer, intent(in) :: nI(nel)
        type(DoubleExc_t), intent(in) :: exc
        logical, intent(in) :: tSign
        HElement_t(dp) :: hel
        unused_var(nI)

        ! Only non-zero contributions if Ms preserved in each term (consider
        ! physical notation).
        hel = sltcnd_2_kernel(exc)

        if (tSign) hel = -hel
    end function sltcnd_2_base

    function sltcnd_2_kernel(exc) result(hel)
        type(DoubleExc_t), intent(in) :: exc
        HElement_t(dp) :: hel
        integer :: id(2, 2), ex(2, 2)
        ! Obtain spatial rather than spin indices if required
        ex = exc%val
        id = gtID(ex)

        if (tReltvy .or. ((G1(ex(1, 1))%Ms == G1(ex(2, 1))%Ms) .and. &
                          (G1(ex(1, 2))%Ms == G1(ex(2, 2))%Ms))) then
            hel = get_umat_el(id(1, 1), id(1, 2), id(2, 1), id(2, 2))
            if (tSpinCorrelator) then
            hel = hel + spinKMatContrib(id(1,1),id(1,2),id(2,1),id(2,2),G1(ex(1,1))%MS,G1(ex(1,2))%MS)
            endif
        else
            hel = (0)
        endif
        if (tReltvy .or. ((G1(ex(1, 1))%Ms == G1(ex(2, 2))%Ms) .and. &
                          (G1(ex(1, 2))%Ms == G1(Ex(2, 1))%Ms))) then
            hel = hel - get_umat_el(id(1, 1), id(1, 2), id(2, 2), id(2, 1))
            if (tSpinCorrelator) then
            hel = hel - spinKMatContrib(id(1,1),id(1,2),id(2,2),id(2,1),G1(ex(1,1))%MS,G1(ex(1,2))%MS)
            endif
        endif

        if (tDampKMat) then
            ! same-spin excitations have only half of KMat acting
            if (G1(ex(1, 1))%Ms == G1(ex(1, 2))%Ms .and. &
  G1(ex(1, 2))%Ms == G1(ex(2, 2))%Ms .and. G1(ex(1, 1))%Ms == G1(ex(2, 1))%Ms) &
     hel = hel + kMatParSpinCorrection(id(1, 1), id(1, 2), id(2, 1), id(2, 2)) &
                 - kMatParSpinCorrection(id(1, 1), id(1, 2), id(2, 2), id(2, 1))
            ! opposite-spin excitations of different orbitals have a different weighting of
            ! normal vs exchange term (not relevant for same-orbtial excits, for lack of exchange)
         if (G1(ex(1, 1))%MS /= G1(ex(1, 2))%MS .and. id(1, 1) /= id(1, 2)) then
                if (G1(ex(1, 1))%MS == G1(ex(2, 1))%MS) &
       hel = hel - kMatOppSpinCorrection(id(1, 1), id(1, 2), id(2, 1), id(2, 2))
                if (G1(ex(1, 1))%MS == G1(ex(2, 2))%MS) &
       hel = hel + kMatOppSpinCorrection(id(1, 1), id(1, 2), id(2, 2), id(2, 1))
            endif
        endif
    end function sltcnd_2_kernel

    !------------------------------------------------------------------------------------------!
    !      slater condon rules for 3-body terms
    !------------------------------------------------------------------------------------------!

    function sltcnd_0_tc(nI, exc) result(hel)
        integer, intent(in) :: nI(nel)
        type(NoExc_t), intent(in) :: exc
        HElement_t(dp) :: hel
        integer :: id(nel)
        integer :: i, j, k

        ! get the diagonal matrix element up to 2nd order
        hel = sltcnd_0_base(nI, exc)
        ! then add the 3-body part
        do i = 1, nel - 2
            do j = i + 1, nel - 1
                do k = j + 1, nel
               hel = hel + get_lmat_el(nI(i), nI(j), nI(k), nI(i), nI(j), nI(k))
                end do
            end do
        end do

    end function sltcnd_0_tc

    function sltcnd_1_tc(nI, ex, tSign) result(hel)
        integer, intent(in) :: nI(nel)
        type(SingleExc_t), intent(in) :: ex
        logical, intent(in) :: tSign
        HElement_t(dp) :: hel
        integer :: i, j

        ! start with the normal matrix element
        hel = sltcnd_1_kernel(nI, ex)

        ! then add the 3-body correction
        do i = 1, nel - 1
            do j = i + 1, nel
                if (ex%val(1) /= nI(i) .and. ex%val(1) /= nI(j)) then
                    hel = hel + get_lmat_el(ex%val(1), nI(i), nI(j), ex%val(2), nI(i), nI(j))
                end if
            end do
        end do

        ! take fermi sign into account
        if (tSign) hel = -hel
    end function sltcnd_1_tc

    function sltcnd_2_tc(nI, exc, tSign) result(hel)
        integer, intent(in) :: nI(nel)
        type(DoubleExc_t), intent(in) :: exc
        logical, intent(in) :: tSign
        HElement_t(dp) :: hel
        integer :: i

        ! get the matrix element up to 2-body terms
        hel = sltcnd_2_kernel(exc)
        ! and the 3-body term
        associate(src1 => exc%val(1, 1), tgt1 => exc%val(2, 1), &
                  src2 => exc%val(1, 2), tgt2 => exc%val(2, 2))
            do i = 1, nel
                if (src1 /= nI(i) .and. src2 /= nI(i)) then
                    hel = hel + get_lmat_el(&
                        src1, src2, nI(i), tgt1, tgt2, nI(i))
                end if
            end do
        end associate

        ! take fermi sign into account
        if (tSign) hel = -hel

    end function sltcnd_2_tc

    function sltcnd_3_tc(ex, tSign) result(hel)
        integer, intent(in) :: ex(2, 3)
        logical, intent(in) :: tSign
        HElement_t(dp) :: hel

        ! this is directly the fully symmetrized entry of the L-matrix
   hel = get_lmat_el(ex(1, 1), ex(1, 2), ex(1, 3), ex(2, 1), ex(2, 2), ex(2, 3))
        ! take fermi sign into account
        if (tSign) hel = -hel
    end function sltcnd_3_tc

    ! dummy function for 3-body matrix elements without tc
    function sltcnd_3_base(ex, tSign) result(hel)
        integer, intent(in) :: ex(2, 3)
        logical, intent(in) :: tSign
        HElement_t(dp) :: hel
        unused_var(ex)
        unused_var(tSign)

        hel = 0
    end function sltcnd_3_base

    !------------------------------------------------------------------------------------------!
    !      slater condon rules for ultracold atoms
    !------------------------------------------------------------------------------------------!

    function sltcnd_0_base_ua(nI, exc) result(hel)
        ! Calculate the  by the SlaterCondon Rules when the two
        ! determinants are the same (so we only need to specify one).
        integer, intent(in) :: nI(nel)
        type(NoExc_t), intent(in) :: exc
        HElement_t(dp) :: hel

        HElement_t(dp) :: hel_sing, hel_doub, hel_tmp
        integer :: id(nel), i, j, idN, idX

        ! Sum in the one electron integrals (KE --> TMAT)
        hel_sing = sum(GetTMATEl(nI, nI))

        ! Obtain the spatial rather than spin indices if required
        id = nI
        !write(6,*) "****",id(:)

        ! Sum in the two electron contributions. Use max(id...) as we cannot
        ! guarantee that if j>i then nI(j)>nI(i).
        hel_doub = (0)
        hel_tmp = (0)
        do i = 1, nel - 1
            do j = i + 1, nel
                idX = max(id(i), id(j))
                idN = min(id(i), id(j))
                hel_doub = hel_doub + get_umat_el(idN, idX, idN, idX)
            enddo
        enddo

        ! Exchange contribution only considered if tExch set.
        ! This is only separated from the above loop to keep "if (tExch)" out
        ! of the tight loop for efficiency.
        if (tExch) then
            do i = 1, nel - 1
                do j = i + 1, nel
                    ! Exchange contribution is zero if I,J are alpha/beta
                    if ((G1(nI(i))%Ms == G1(nI(j))%Ms) .or. tReltvy) then
                        idX = max(id(i), id(j))
                        idN = min(id(i), id(j))
                        hel_tmp = hel_tmp - get_umat_el(idN, idX, idX, idN)
                    endif
                enddo
            enddo
        endif
        hel = hel_doub + hel_tmp + hel_sing

    end function sltcnd_0_base_ua

    function sltcnd_1_base_ua(nI, ex, tSign) result(hel)
        ! Calculate the  by the Slater-Condon Rules when the two
        ! determinants differ by one orbital exactly.
        integer, intent(in) :: nI(nel)
        type(SingleExc_t), intent(in) :: ex
        logical, intent(in) :: tSign
        HElement_t(dp) :: hel

        ! Sum in the diagonal terms (same in both dets)
        ! Coulomb term only included if Ms values of ex(1) and ex(2) are the
        ! same.

        hel = sltcnd_1_kernel_ua(nI, ex)
        if (tSign) hel = -hel
    end function sltcnd_1_base_ua

    function sltcnd_1_kernel_ua(nI, exc) result(hel)
        integer, intent(in) :: nI(nel)
        type(SingleExc_t) :: exc
        HElement_t(dp) :: hel
        integer :: i, id, id_ex(2)

        ! Obtain spatial rather than spin indices if required
        id_ex = exc%val

        hel = (0)
        if (tReltvy .or. (G1(exc%val(1))%Ms == G1(exc%val(2))%Ms)) then
            do i = 1, nel
                if (exc%val(1) /= nI(i)) then
                    id = nI(i)
                    hel = hel + get_umat_el(id_ex(1), id, id_ex(2), id)
                endif
            enddo
        endif
        ! Exchange contribution is only considered if tExch set.
        ! This is only separated from the above loop to keep "if (tExch)" out
        ! of the tight loop for efficiency.
        if (tExch .and. ((G1(exc%val(1))%Ms == G1(exc%val(2))%Ms) .or. tReltvy)) then
            do i = 1, nel
                if (exc%val(1) /= nI(i)) then
                    if (tReltvy .or. (G1(exc%val(1))%Ms == G1(nI(i))%Ms)) then
                        id = nI(i)
                        hel = hel - get_umat_el(id_ex(1), id, id, id_ex(2))
                    endif
                endif
            enddo
        endif
        ! consider the non-diagonal part of the kinetic energy -
        ! <psi_a|T|psi_a'> where a, a' are the only basis fns that differ in
        ! nI, nJ
        hel = hel + GetTMATEl(exc%val(1), exc%val(2))
    end function sltcnd_1_kernel_ua

    function sltcnd_2_base_ua(nI, ex, tSign) result(hel)
        ! Calculate the  by the Slater-Condon Rules when the two
        ! determinants differ by two orbitals exactly (the simplest case).
        integer, intent(in) :: nI(nel)
        type(DoubleExc_t), intent(in) :: ex
        logical, intent(in) :: tSign
        HElement_t(dp) :: hel
        unused_var(nI)

        ! Only non-zero contributions if Ms preserved in each term (consider
        ! physical notation).
        hel = sltcnd_2_kernel_ua(ex)

        if (tSign) hel = -hel
    end function sltcnd_2_base_ua

    function sltcnd_2_kernel_ua(ex) result(hel)
        implicit none
        type(DoubleExc_t), intent(in) :: ex
        HElement_t(dp) :: hel
        integer :: id(2, 2)
        ! Obtain spatial rather than spin indices if required
        id = ex%val

        associate(src1 => ex%val(1, 1), tgt1 => ex%val(2, 1), &
                  src2 => ex%val(1, 2), tgt2 => ex%val(2, 2))

            if (tReltvy .or. ((G1(src1)%Ms == G1(tgt1)%Ms) .and. &
                              (G1(src2)%Ms == G1(tgt2)%Ms))) then
                hel = get_umat_el(id(1, 1), id(1, 2), id(2, 1), id(2, 2))
            else
                hel = (0)
            endif
            if (tReltvy .or. ((G1(src1)%Ms == G1(tgt2)%Ms) .and. &
                              (G1(src2)%Ms == G1(tgt1)%Ms))) then
                hel = hel - get_umat_el(id(1, 1), id(1, 2), id(2, 2), id(2, 1))
            endif
        end associate
    end function sltcnd_2_kernel_ua

    function sltcnd_2_kernel_ua_3b(nI, exc) result(hel)
        implicit none
        integer, intent(in) :: nI(nel)
        type(DoubleExc_t), intent(in) :: exc
        HElement_t(dp) :: hel
        integer :: id(2, 2), ex(2, 2)
        ! Obtain spatial rather than spin indices if required
        id = exc%val
        ex = exc%val

        hel = (0)
        if (G1(ex(1, 1))%Ms == G1(ex(1, 2))%Ms) then
            if (tReltvy .or. ((G1(ex(2, 1))%Ms == G1(ex(2, 2))%Ms) .and. &
                              (G1(ex(1, 1))%Ms == G1(ex(2, 2))%Ms))) then
     hel = get_contact_umat_el_3b_sp(id(1, 1), id(1, 2), id(2, 1), id(2, 2)) - &
               get_contact_umat_el_3b_sp(id(1, 1), id(1, 2), id(2, 2), id(2, 1))
            endif
        else
            ! We have an additional sign factor due to the exchange of the creation
            ! operators:
            !a_(p-k)^+ a_(s+k)^+ a_q^+ a_q a_s a_p -> -a_q^+ a_(s+k)^+ a_(p-k)^+ a_q a_s a_p
            !a_(p-k)^+ a_(s+k)^+ a_q^+ a_q a_s a_p -> -a_(p-k)^+ a_q^+ a_(s+k)^+ a_q a_s a_p
            if (tReltvy .or. ((G1(ex(1, 1))%Ms == G1(ex(2, 1))%Ms) .and. &
                              (G1(ex(1, 2))%Ms == G1(ex(2, 2))%Ms))) then
   hel = -get_contact_umat_el_3b_sap(id(1, 1), id(1, 2), id(2, 1), id(2, 2), nI)
            endif
            if (tReltvy .or. ((G1(ex(1, 1))%Ms == G1(ex(2, 2))%Ms) .and. &
                              (G1(ex(1, 2))%Ms == G1(Ex(2, 1))%Ms))) then
           hel = hel + get_contact_umat_el_3b_sap (id(1,1), id(1,2), id(2,2), id(2,1), nI)
            endif
        endif

    end function sltcnd_2_kernel_ua_3b

    function sltcnd_0_tc_ua(nI, exc) result(hel)
        integer, intent(in) :: nI(nel)
        type(NoExc_t), intent(in) :: exc
        HElement_t(dp) :: hel
        integer :: id(nel)
        integer :: i, j, k

        ! get the diagonal matrix element up to 2nd order
        hel = sltcnd_0_base_ua(nI, exc)
        ! then add the 3-body part
        do i = 1, nel - 2
            do j = i + 1, nel - 1
                do k = j + 1, nel
            hel = hel + get_lmat_el_ua(nI(i), nI(j), nI(k), nI(i), nI(j), nI(k))
                end do
            end do
        end do

    end function sltcnd_0_tc_ua

    function sltcnd_1_tc_ua(nI, exc, tSign) result(hel)
        integer, intent(in) :: nI(nel)
        type(SingleExc_t), intent(in) :: exc
        logical, intent(in) :: tSign
        HElement_t(dp) :: hel
        integer :: i, j

        ! start with the normal matrix element
        hel = sltcnd_1_kernel_ua(nI, exc)

        ! then add the 3-body correction
        do i = 1, nel - 1
            do j = i + 1, nel
                if (exc%val(1) /= nI(i) .and. exc%val(1) /= nI(j)) then
            hel = hel + get_lmat_el_ua(exc%val(1), nI(i), nI(j), exc%val(2), nI(i), nI(j))
!      print *, "from", ex(1), nI(i),nI(j)
!      print *, "to", ex(2),nI(i),nI(j)
!      print *, "hel", hel,get_lmat_el_ua(ex(1),nI(i),nI(j),ex(2),nI(i),nI(j))
                endif
            end do
        end do

        ! take fermi sign into account
        if (tSign) hel = -hel
    end function sltcnd_1_tc_ua

    function sltcnd_2_tc_ua(nI, ex, tSign) result(hel)
        integer, intent(in) :: nI(nel)
        type(DoubleExc_t), intent(in) :: ex
        logical, intent(in) :: tSign
        HElement_t(dp) :: hel, heltc
        integer :: i

        ! get the matrix element up to 2-body terms
        hel = sltcnd_2_kernel_ua(ex)
        ! and the 3-body term
!       heltc= 0.d0
!     do i = 1, nel
!        if(ex(1,1).ne.nI(i) .and. ex(1,2).ne.nI(i)) then! &
!        heltc = heltc + get_lmat_el_ua(ex(1,1),ex(1,2),nI(i),ex(2,1),ex(2,2),nI(i))
!       endif
!     end do
!         hel=heltc
!         heltc=0.d0
        heltc = sltcnd_2_kernel_ua_3b(nI, ex)
!       if(dabs(hel-heltc).gt.0.000000001) then
!       write(6,*) 'nI', nI(1:nel)
!       write(6,*) 'ex', ex(1,1:2), '->', ex(2,1:2)
!       write(6,*) 'heltc', heltc
!       write(6,*) 'hel', hel
!               call stop_all()
!       endif

        hel = hel + heltc

        ! take fermi sign into account
        if (tSign) hel = -hel

    end function sltcnd_2_tc_ua

    function sltcnd_3_tc_ua(ex, tSign) result(hel)
        integer, intent(in) :: ex(2, 3)
        logical, intent(in) :: tSign
        HElement_t(dp) :: hel

        ! this is directly the fully symmetrized entry of the L-matrix
hel = get_lmat_el_ua(ex(1, 1), ex(1, 2), ex(1, 3), ex(2, 1), ex(2, 2), ex(2, 3))
        ! take fermi sign into account
        if (tSign) hel = -hel
    end function sltcnd_3_tc_ua

    HElement_t(dp) function sltcnd_excit_NoExc_t(ref, exc)
        integer, intent(in) :: ref(nel)
        type(NoExc_t), intent(in), optional :: exc

        sltcnd_excit_NoExc_t = sltcnd_0(ref, exc)
    end function

    HElement_t(dp) function sltcnd_excit_SingleExc_t(ref, exc, tParity)
        integer, intent(in) :: ref(nel)
        type(SingleExc_t), intent(in) :: exc
        logical, intent(in) :: tParity

        sltcnd_excit_SingleExc_t = sltcnd_1(ref, exc, tParity)
    end function

    HElement_t(dp) function sltcnd_excit_DoubleExc_t(ref, exc, tParity)
        integer, intent(in) :: ref(nel)
        type(DoubleExc_t), intent(in) :: exc
        logical, intent(in) :: tParity

        sltcnd_excit_DoubleExc_t = sltcnd_2(ref, exc, tParity)
    end function

    HElement_t(dp) function sltcnd_excit_TripleExc_t(ref, exc, tParity)
        integer, intent(in) :: ref(nel)
        type(TripleExc_t), intent(in) :: exc
        logical, intent(in) :: tParity
        @:unused_var(ref)

        sltcnd_excit_TripleExc_t = sltcnd_3(exc%val, tParity)
    end function

end module

