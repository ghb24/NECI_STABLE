#include "macros.h"
#:include "macros.fpph"
#:set max_excit_rank = 3
#:set excit_ranks = list(range(max_excit_rank + 1))
#:set excitations = [f'Excite_{i}_t' for i in ['Further'] + excit_ranks]
#:set trivial_excitations = excitations[:2]
#:set non_trivial_excitations = excitations[2:]
#:set classic_abinit_excitations = excitations[1:max_excit_rank + 1]

!>  @brief
!>      A module to evaluate the Slater-Condon Rules.
!>
!>  @details
!>  Heavily relies on the excitation_types module to represent different excitations.
!>
!>  The main functions are sltcnd_excit and dyn_sltcnd_excit.
!>  Because sltcnd_excit dispatches statically at compile time,
!>  depending on the excitation type, it is the preferred function,
!>  if the excitation level is known at compile time.
!>
!>  If the excitation level is not known at compile time,
!>  use dyn_sltcnd_excit which accepts a polymorphic excitation_t class.
!>
!>  The procedures create_excitation, get_excitation, and get_bit_excitation
!>  from the excitation_types module
!>  can be used, to create excitations from nIs, or iluts at runtime.
module sltcnd_mod
    use SystemData, only: nel, nBasisMax, tExch, G1, ALAT, tReltvy, &
                          t_mol_3_body, t_ueg_3_body, tContact
    use SystemData, only: nBasis!, iSpinSkip
    ! HACK - We use nBasisMax(2,3) here rather than iSpinSkip, as it appears
    !        to be more reliably set (see for example test H2O_RI)
    ! TODO: We need to sort this out so that they are consistent
    !       --> Talk to George/Alex to see what impact that might have?
    use constants, only: dp, n_int, maxExcit
    use UMatCache, only: GTID, UMatInd
    use IntegralsData, only: UMAT, t_use_tchint_lib
    use OneEInts, only: GetTMatEl, TMat2D
    use procedure_pointers, only: get_umat_el
    use excitation_types, only: excitation_t, Excite_0_t, Excite_1_t, Excite_2_t, &
                                Excite_3_t, UNKNOWN, get_excitation, get_bit_excitation, &
                                create_excitation, Excite_Further_t
    use orb_idx_mod, only: SpinOrbIdx_t
    use DetBitOps, only: count_open_orbs, FindBitExcitLevel
    use bit_reps, only: NIfTot
    use LMat_mod, only: get_lmat_el, get_lmat_el_ua, external_lMat_matel
    use gen_coul_ueg_mod, only: get_contact_umat_el_3b_sp, get_contact_umat_el_3b_sap
    use SD_spin_purification_mod, only: possible_purification_methods, SD_spin_purification, &
                spin_pure_J, S2_expval_exc, ladder_op_exc
    use util_mod, only: stop_all

    implicit none
    private
    public :: initSltCndPtr, &
              ! statically known
              sltcnd_excit, sltcnd_2_kernel, &
              ! dynamically known
              dyn_sltcnd_excit_old, dyn_sltcnd_excit, &
              sltcnd_compat, sltcnd, sltcnd_knowIC, &
              CalcFockOrbEnergy, sumfock, sltcnd_0_base, sltcnd_0_tc

!>  @brief
!>      Evaluate Matrix Element for different excitations
!>      using the Slater-Condon rules.
!>
!>  @details
!>  This generic function uses compile time dispatch.
!>  This means that exc cannot be just of class(excitation_t) but has
!>  to be a proper non-polymorphic type.
!>  For run time dispatch use dyn_sltcnd_excit.
!>
!>  @param[in] exc, An excitation of a subtype of excitation_t.
    interface sltcnd_excit
    #:for excitation_t in excitations
        module procedure sltcnd_excit_${excitation_t}$
    #:endfor
        module procedure sltcnd_excit_SpinOrbIdx_t_Excite_1_t
        module procedure sltcnd_excit_SpinOrbIdx_t_Excite_2_t
    end interface

    abstract interface
        #:for rank in excit_ranks
        HElement_t(dp) function sltcnd_${rank}$_t(nI, exc, tParity) result(hel)
            import :: dp, nel, Excite_${rank}$_t
            integer, intent(in) :: nI(nel)
            type(Excite_${rank}$_t), intent(in) :: exc
            logical, intent(in) :: tParity
        end function sltcnd_${rank}$_t
        #:endfor
    end interface

    #:for rank in excit_ranks
    procedure(sltcnd_${rank}$_t), pointer :: sltcnd_${rank}$ => null()
    #:endfor

contains

    subroutine initSltCndPtr()
        use SystemData, only: tSmallBasisForThreeBody
        character(*), parameter :: this_routine = 'initSltCndPtr'

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
            if (t_mol_3_body .or. t_ueg_3_body .and. nel >= 2) then
                sltcnd_0 => sltcnd_0_tc
                sltcnd_1 => sltcnd_1_tc
                sltcnd_2 => sltcnd_2_tc
                sltcnd_3 => sltcnd_3_tc
            else if (allocated(SD_spin_purification)) then
                if (SD_spin_purification == possible_purification_methods%TRUNCATED_LADDER) then
                    sltcnd_0 => sltcnd_0_base
                else if (SD_spin_purification == possible_purification_methods%ONLY_LADDER) then
                    sltcnd_0 => sltcnd_0_purify_spin_only_ladder
                else if (SD_spin_purification == possible_purification_methods%FULL_S2) then
                    sltcnd_0 => sltcnd_0_purify_spin_full_s2
                else
                    call stop_all(this_routine, 'Invalid options for SD_spin_purification')
                end if
                sltcnd_2 => sltcnd_2_purify_spin

                ! Unaffected by < I | S^2 | J >
                sltcnd_1 => sltcnd_1_base
                sltcnd_3 => sltcnd_3_base

            else
                sltcnd_0 => sltcnd_0_base
                sltcnd_1 => sltcnd_1_base
                sltcnd_2 => sltcnd_2_base
                sltcnd_3 => sltcnd_3_base
            end if

        end if
    end subroutine initSltCndPtr

!>  @brief
!>      Evaluate Matrix Element for different excitations
!>      using the Slater-Condon rules.
!>
!>  @details
!>  This generic function uses run time dispatch.
!>  This means that exc can be any subtype of class(excitation_t).
!>  For performance reason it is advised to use sltcnd_excit,
!>  if the actual type is known at compile time.
!>
!>  @param[in] ref, The reference determinant as array of occupied orbital indices.
!>  @param[in] exc, An excitation of type excitation_t.
!>  @param[in] tParity, The parity of the excitation.
    function dyn_sltcnd_excit(ref, exc, tParity) result(hel)
        integer, intent(in) :: ref(nel)
        class(excitation_t), intent(in) :: exc
        logical, intent(in) :: tParity
        HElement_t(dp) :: hel
        character(*), parameter :: this_routine = 'dyn_sltcnd_excit'

        ! The compiler has to statically know, of what type exc is.
        select type (exc)
        #:for Excitation_t in excitations
        type is (${Excitation_t}$)
            hel = sltcnd_excit(ref, exc, tParity)
        #:endfor
        class default
            call stop_all(this_routine, "Error in downcast.")
        end select
    end function

    function dyn_sltcnd_excit_old(nI, IC, ex, tParity) result(hel)

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

        class(excitation_t), allocatable :: exc

        if (IC /= 0 .and. .not. (present(ex) .and. present(tParity))) &
            call stop_all(this_routine, "ex and tParity must be provided to &
                          &sltcnd_excit for all IC /= 0")
        call create_excitation(exc, IC, ex)
        hel = dyn_sltcnd_excit(nI, exc, tParity)
    end function

    function sltcnd_compat(nI, nJ, IC) result(hel)
        integer, intent(in) :: nI(nel), nJ(nel), IC
        HElement_t(dp) :: hel

        class(excitation_t), allocatable :: exc
        logical :: tParity

        call get_excitation(nI, nJ, IC, exc, tParity)
        hel = dyn_sltcnd_excit(nI, exc, tParity)
    end function sltcnd_compat

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

        class(excitation_t), allocatable :: exc
        logical :: tParity

        call get_bit_excitation(ilutI, ilutJ, IC, exc, tParity)
        hel = dyn_sltcnd_excit(nI, exc, tParity)

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
        end do

        ! Exchange contribution only considered if tExch set.
        ! This is only separated from the above loop to keep "if (tExch)" out
        ! of the tight loop for efficiency.
        do j = 1, nel
            ! Exchange contribution is zero if I,J are alpha/beta
            if (tReltvy .or. (G1(Orb)%Ms == G1(HFDet(j))%Ms)) then
                idN = idHF(j)
                hel = hel - get_umat_el(idOrb, idN, idN, idOrb)
            end if
        end do
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
            end do
        end do

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
                    end if
                end do
            end do
        end if
        hel = hel_doub + hel_tmp + hel_sing

    end function SumFock

    pure function sltcnd_0_base(nI, exc, tSign) result(hel)
        ! Calculate the  by the SlaterCondon Rules when the two
        ! determinants are the same (so we only need to specify one).
        integer, intent(in) :: nI(nel)
        type(Excite_0_t), intent(in) :: exc
        HElement_t(dp) :: hel, hel_sing, hel_doub, hel_tmp
        integer :: id(nel), i, j
        logical, intent(in) :: tSign

        @:unused_var(exc, tSign)

        ! Sum in the one electron integrals (KE --> TMAT)
        hel_sing = sum(GetTMATEl(nI, nI))

        ! Obtain the spatial rather than spin indices if required
        id = gtID(nI)

        ! Sum in the two electron contributions.
        hel_doub = (0)
        hel_tmp = (0)
        do i = 1, nel - 1
            do j = i + 1, nel
                hel_doub = hel_doub + get_umat_el(id(i), id(j), id(i), id(j))
            end do
        end do

        ! Exchange contribution only considered if tExch set.
        ! This is only separated from the above loop to keep "if (tExch)" out
        ! of the tight loop for efficiency.
        if (tExch) then
            do i = 1, nel - 1
                do j = i + 1, nel
                    ! Exchange contribution is zero if I,J are alpha/beta
                    if ((G1(nI(i))%Ms == G1(nI(j))%Ms) .or. tReltvy) then
                        hel_tmp = hel_tmp - get_umat_el(id(i), id(j), id(j), id(i))
                    end if
                end do
            end do
        end if
        hel = hel_doub + hel_tmp + hel_sing

    end function sltcnd_0_base

    function sltcnd_1_base(nI, ex, tSign) result(hel)
        ! Calculate the  by the Slater-Condon Rules when the two
        ! determinants differ by one orbital exactly.
        integer, intent(in) :: nI(nel)
        type(Excite_1_t), intent(in) :: ex
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
        type(Excite_1_t), intent(in) :: ex
        HElement_t(dp) :: hel
        integer :: i, id, id_ex(2)

        ! Obtain spatial rather than spin indices if required
        id_ex = gtID(ex%val(:,1))

        hel = (0)
        if (tReltvy .or. (G1(ex%val(1, 1))%Ms == G1(ex%val(2, 1))%Ms)) then
            do i = 1, nel
                if (ex%val(1, 1) /= nI(i)) then
                    id = gtID(nI(i))
                    hel = hel + get_umat_el(id_ex(1), id, id_ex(2), id)
                end if
            end do
        end if
        ! Exchange contribution is only considered if tExch set.
        ! This is only separated from the above loop to keep "if (tExch)" out
        ! of the tight loop for efficiency.
        if (tExch .and. ((G1(ex%val(1, 1))%Ms == G1(ex%val(2, 1))%Ms) .or. tReltvy)) then
            do i = 1, nel
                if (ex%val(1, 1) /= nI(i)) then
                    if (tReltvy .or. (G1(ex%val(1, 1))%Ms == G1(nI(i))%Ms)) then
                        id = gtID(nI(i))
                        hel = hel - get_umat_el(id_ex(1), id, id, id_ex(2))
                    end if
                end if
            end do
        end if
        ! consider the non-diagonal part of the kinetic energy -
        ! <psi_a|T|psi_a'> where a, a' are the only basis fns that differ in
        ! nI, nJ
        hel = hel + GetTMATEl(ex%val(1, 1), ex%val(2, 1))
    end function sltcnd_1_kernel

    function sltcnd_2_base(nI, exc, tSign) result(hel)
        ! Calculate the  by the Slater-Condon Rules when the two
        ! determinants differ by two orbitals exactly (the simplest case).
        integer, intent(in) :: nI(nel)
        type(Excite_2_t), intent(in) :: exc
        logical, intent(in) :: tSign
        HElement_t(dp) :: hel
        @:unused_var(nI)

        ! Only non-zero contributions if Ms preserved in each term (consider
        ! physical notation).
        hel = sltcnd_2_kernel(exc)

        if (tSign) hel = -hel
    end function sltcnd_2_base

    function sltcnd_2_kernel(exc) result(hel)
        type(Excite_2_t), intent(in) :: exc
        HElement_t(dp) :: hel
        integer :: id(2, 2), ex(2, 2)
        ! Obtain spatial rather than spin indices if required
        ex = exc%val
        id = gtID(ex)

        if (tReltvy .or. ((G1(ex(1, 1))%Ms == G1(ex(2, 1))%Ms) .and. &
                          (G1(ex(1, 2))%Ms == G1(ex(2, 2))%Ms))) then
            hel = get_umat_el(id(1, 1), id(1, 2), id(2, 1), id(2, 2))
        else
            hel = (0)
        end if
        if (tReltvy .or. ((G1(ex(1, 1))%Ms == G1(ex(2, 2))%Ms) .and. &
                          (G1(ex(1, 2))%Ms == G1(Ex(2, 1))%Ms))) then
            hel = hel - get_umat_el(id(1, 1), id(1, 2), id(2, 2), id(2, 1))
        end if

    end function sltcnd_2_kernel

    !------------------------------------------------------------------------------------------!
    !      slater condon rules for 3-body terms
    !------------------------------------------------------------------------------------------!

    function sltcnd_0_tc(nI, exc, tParity) result(hel)
        integer, intent(in) :: nI(nel)
        type(Excite_0_t), intent(in) :: exc
        HElement_t(dp) :: hel
        integer :: i, j, k
        integer :: dummy(1,0)
        logical, intent(in) :: tParity

        ! get the diagonal matrix element up to 2nd order
        hel = sltcnd_0_base(nI, exc, tParity)
        ! then add the 3-body part
        if(t_use_tchint_lib) then
            hel = hel + external_lMat_matel(nI, dummy)
        else
            do i = 1, nel - 2
                do j = i + 1, nel - 1
                    do k = j + 1, nel
                        hel = hel + get_lmat_el(nI(i), nI(j), nI(k), nI(i), nI(j), nI(k))
                    end do
                end do
            end do
        end if
    end function sltcnd_0_tc

    function sltcnd_1_tc(nI, ex, tSign) result(hel)
        integer, intent(in) :: nI(nel)
        type(Excite_1_t), intent(in) :: ex
        logical, intent(in) :: tSign
        HElement_t(dp) :: hel
        integer :: i, j

        ! start with the normal matrix element
        hel = sltcnd_1_kernel(nI, ex)

        ! then add the 3-body correction
        if(t_use_tchint_lib) then
            hel = hel + external_lMat_matel(nI, reshape(ex%val,(/2,1/)))
        else
            do i = 1, nel - 1
                do j = i + 1, nel
                    if (ex%val(1, 1) /= nI(i) .and. ex%val(1, 1) /= nI(j)) then
                        hel = hel + get_lmat_el(ex%val(1, 1), nI(i), nI(j), ex%val(2, 1), nI(i), nI(j))
                    end if
                end do
            end do
        end if
        ! take fermi sign into account
        if (tSign) hel = -hel
    end function sltcnd_1_tc

    function sltcnd_2_tc(nI, exc, tSign) result(hel)
        integer, intent(in) :: nI(nel)
        type(Excite_2_t), intent(in) :: exc
        logical, intent(in) :: tSign
        HElement_t(dp) :: hel
        integer :: i

        ! get the matrix element up to 2-body terms
        hel = sltcnd_2_kernel(exc)

        ! and the 3-body term
        if(t_use_tchint_lib) then
            hel = hel + external_lMat_matel(nI, exc%val)
        else
            associate(src1 => exc%val(1, 1), tgt1 => exc%val(2, 1), &
                src2 => exc%val(1, 2), tgt2 => exc%val(2, 2))
                do i = 1, nel
                    if (src1 /= nI(i) .and. src2 /= nI(i)) then
                    hel = hel + get_lmat_el( &
                        src1, src2, nI(i), tgt1, tgt2, nI(i))
                    end if
                end do
            end associate
        end if
        ! take fermi sign into account
        if (tSign) hel = -hel

    end function sltcnd_2_tc

    function sltcnd_3_tc(nI, ex, tSign) result(hel)
        type(Excite_3_t), intent(in) :: ex
        logical, intent(in) :: tSign
        HElement_t(dp) :: hel
        integer :: dummy(1)
        integer, intent(in) :: nI(nel)
        @:unused_var(nI)

        ! this is directly the fully symmetrized entry of the L-matrix
        if(t_use_tchint_lib) then
            hel = external_lMat_matel(dummy, ex%val)
        else
            associate(ex => ex%val)
                hel = get_lmat_el(ex(1, 1), ex(1, 2), ex(1, 3), &
                                  ex(2, 1), ex(2, 2), ex(2, 3))
            end associate
        endif
        ! take fermi sign into account
        if (tSign) hel = -hel
    end function sltcnd_3_tc

    ! dummy function for 3-body matrix elements without tc
    function sltcnd_3_base(nI, ex, tSign) result(hel)
        integer, intent(in) :: nI(nel)
        type(Excite_3_t), intent(in) :: ex
        logical, intent(in) :: tSign
        HElement_t(dp) :: hel
        @:unused_var(ex, tSign, nI)
        hel = 0
    end function sltcnd_3_base

    !------------------------------------------------------------------------------------------!
    !      slater condon rules for ultracold atoms
    !------------------------------------------------------------------------------------------!

    pure function sltcnd_0_base_ua(nI, exc, tSign) result(hel)
        ! Calculate the  by the SlaterCondon Rules when the two
        ! determinants are the same (so we only need to specify one).
        integer, intent(in) :: nI(nel)
        type(Excite_0_t), intent(in) :: exc
        logical, intent(in) :: tSign
        HElement_t(dp) :: hel

        HElement_t(dp) :: hel_sing, hel_doub, hel_tmp
        integer :: id(nel), i, j, idN, idX

        @:unused_var(exc, tSign)

        ! Sum in the one electron integrals (KE --> TMAT)
        hel_sing = sum(GetTMATEl(nI, nI))

        ! Obtain the spatial rather than spin indices if required
        id = nI

        ! Sum in the two electron contributions. Use max(id...) as we cannot
        ! guarantee that if j>i then nI(j)>nI(i).
        hel_doub = (0)
        hel_tmp = (0)
        do i = 1, nel - 1
            do j = i + 1, nel
                idX = max(id(i), id(j))
                idN = min(id(i), id(j))
                hel_doub = hel_doub + get_umat_el(idN, idX, idN, idX)
            end do
        end do

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
                    end if
                end do
            end do
        end if
        hel = hel_doub + hel_tmp + hel_sing

    end function sltcnd_0_base_ua

    function sltcnd_1_base_ua(nI, ex, tSign) result(hel)
        ! Calculate the  by the Slater-Condon Rules when the two
        ! determinants differ by one orbital exactly.
        integer, intent(in) :: nI(nel)
        type(Excite_1_t), intent(in) :: ex
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
        type(Excite_1_t) :: exc
        HElement_t(dp) :: hel
        integer :: i, id, id_ex(2)

        ! Obtain spatial rather than spin indices if required
        id_ex = exc%val(:, 1)

        hel = (0)
        if (tReltvy .or. (G1(exc%val(1, 1))%Ms == G1(exc%val(2, 1))%Ms)) then
            do i = 1, nel
                if (exc%val(1, 1) /= nI(i)) then
                    id = nI(i)
                    hel = hel + get_umat_el(id_ex(1), id, id_ex(2), id)
                end if
            end do
        end if
        ! Exchange contribution is only considered if tExch set.
        ! This is only separated from the above loop to keep "if (tExch)" out
        ! of the tight loop for efficiency.
        if (tExch .and. ((G1(exc%val(1, 1))%Ms == G1(exc%val(2, 1))%Ms) .or. tReltvy)) then
            do i = 1, nel
                if (exc%val(1, 1) /= nI(i)) then
                    if (tReltvy .or. (G1(exc%val(1, 1))%Ms == G1(nI(i))%Ms)) then
                        id = nI(i)
                        hel = hel - get_umat_el(id_ex(1), id, id, id_ex(2))
                    end if
                end if
            end do
        end if
        ! consider the non-diagonal part of the kinetic energy -
        ! <psi_a|T|psi_a'> where a, a' are the only basis fns that differ in
        ! nI, nJ
        hel = hel + GetTMATEl(exc%val(1, 1), exc%val(2, 1))
    end function sltcnd_1_kernel_ua

    function sltcnd_2_base_ua(nI, ex, tSign) result(hel)
        ! Calculate the  by the Slater-Condon Rules when the two
        ! determinants differ by two orbitals exactly (the simplest case).
        integer, intent(in) :: nI(nel)
        type(Excite_2_t), intent(in) :: ex
        logical, intent(in) :: tSign
        HElement_t(dp) :: hel

        @:unused_var(nI)

        ! Only non-zero contributions if Ms preserved in each term (consider
        ! physical notation).
        hel = sltcnd_2_kernel_ua(ex)

        if (tSign) hel = -hel
    end function sltcnd_2_base_ua

    function sltcnd_2_kernel_ua(ex) result(hel)
        implicit none
        type(Excite_2_t), intent(in) :: ex
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
            end if
            if (tReltvy .or. ((G1(src1)%Ms == G1(tgt2)%Ms) .and. &
                              (G1(src2)%Ms == G1(tgt1)%Ms))) then
                hel = hel - get_umat_el(id(1, 1), id(1, 2), id(2, 2), id(2, 1))
            end if
        end associate
    end function sltcnd_2_kernel_ua

    function sltcnd_2_kernel_ua_3b(nI, exc) result(hel)
        implicit none
        integer, intent(in) :: nI(nel)
        type(Excite_2_t), intent(in) :: exc
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
            end if
        else
            ! We have an additional sign factor due to the exchange of the creation
            ! operators:
            !a_(p-k)^+ a_(s+k)^+ a_q^+ a_q a_s a_p -> -a_q^+ a_(s+k)^+ a_(p-k)^+ a_q a_s a_p
            !a_(p-k)^+ a_(s+k)^+ a_q^+ a_q a_s a_p -> -a_(p-k)^+ a_q^+ a_(s+k)^+ a_q a_s a_p
            if (tReltvy .or. ((G1(ex(1, 1))%Ms == G1(ex(2, 1))%Ms) .and. &
                              (G1(ex(1, 2))%Ms == G1(ex(2, 2))%Ms))) then
                hel = -get_contact_umat_el_3b_sap(id(1, 1), id(1, 2), id(2, 1), id(2, 2), nI)
            end if
            if (tReltvy .or. ((G1(ex(1, 1))%Ms == G1(ex(2, 2))%Ms) .and. &
                              (G1(ex(1, 2))%Ms == G1(Ex(2, 1))%Ms))) then
                hel = hel + get_contact_umat_el_3b_sap(id(1, 1), id(1, 2), id(2, 2), id(2, 1), nI)
            end if
        end if

    end function sltcnd_2_kernel_ua_3b

    pure function sltcnd_0_tc_ua(nI, exc, tParity) result(hel)
        integer, intent(in) :: nI(nel)
        type(Excite_0_t), intent(in) :: exc
        HElement_t(dp) :: hel
        integer :: i, j, k
        logical, intent(in) :: tParity

        @:unused_var(tParity)

        ! get the diagonal matrix element up to 2nd order
        hel = sltcnd_0_base_ua(nI, exc, tParity)
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
        type(Excite_1_t), intent(in) :: exc
        logical, intent(in) :: tSign
        HElement_t(dp) :: hel
        integer :: i, j

        ! start with the normal matrix element
        hel = sltcnd_1_kernel_ua(nI, exc)

        ! then add the 3-body correction
        do i = 1, nel - 1
            do j = i + 1, nel
                if (exc%val(1, 1) /= nI(i) .and. exc%val(1, 1) /= nI(j)) then
                    hel = hel + get_lmat_el_ua(exc%val(1, 1), nI(i), nI(j), exc%val(2, 1), nI(i), nI(j))
                end if
            end do
        end do

        ! take fermi sign into account
        if (tSign) hel = -hel
    end function sltcnd_1_tc_ua

    function sltcnd_2_tc_ua(nI, ex, tSign) result(hel)
        integer, intent(in) :: nI(nel)
        type(Excite_2_t), intent(in) :: ex
        logical, intent(in) :: tSign
        HElement_t(dp) :: hel, heltc

        ! get the matrix element up to 2-body terms
        hel = sltcnd_2_kernel_ua(ex)
        ! and the 3-body term
        heltc = sltcnd_2_kernel_ua_3b(nI, ex)
        hel = hel + heltc

        ! take fermi sign into account
        if (tSign) hel = -hel
    end function sltcnd_2_tc_ua

    function sltcnd_3_tc_ua(nI, ex, tSign) result(hel)
        type(Excite_3_t), intent(in) :: ex
        logical, intent(in) :: tSign
        HElement_t(dp) :: hel
        integer, intent(in) :: nI(nEl)

        @:unused_var(nI)

        ! this is directly the fully symmetrized entry of the L-matrix
        associate(ex => ex%val)
            hel = get_lmat_el_ua(ex(1, 1), ex(1, 2), ex(1, 3), &
                                 ex(2, 1), ex(2, 2), ex(2, 3))
        end associate
        ! take fermi sign into account
        if (tSign) hel = -hel
    end function sltcnd_3_tc_ua

    function sltcnd_0_purify_spin_only_ladder(nI, exc, tParity) result(hel)
        integer, intent(in) :: nI(nel)
        type(Excite_0_t), intent(in) :: exc
        logical, intent(in) :: tParity
        HElement_t(dp) :: hel
        hel = sltcnd_0_base(nI, exc, tParity) + spin_pure_J * ladder_op_exc(nI, exc)
    end function

    function sltcnd_0_purify_spin_full_s2(nI, exc, tParity) result(hel)
        integer, intent(in) :: nI(nel)
        type(Excite_0_t), intent(in) :: exc
        logical, intent(in) :: tParity
        HElement_t(dp) :: hel
        hel = sltcnd_0_base(nI, exc, tParity) + spin_pure_J * S2_expval_exc(nI, exc)
    end function


    function sltcnd_2_purify_spin(nI, exc, tSign) result(hel)
        integer, intent(in) :: nI(nel)
        type(Excite_2_t), intent(in) :: exc
        logical, intent(in) :: tSign
        HElement_t(dp) :: hel
        hel = sltcnd_2_base(nI, exc, tSign) + spin_pure_J * S2_expval_exc(nI, exc)
    end function

    !>  @brief
    !>      Evaluate Matrix Element for the reference (no Excitation).
    !>
    !>  @param[in] ref, An array of occupied orbital indices in the reference.
    !>  @param[in] exc, An excitation of type Excite_0_t.
    HElement_t(dp) function sltcnd_excit_Excite_0_t(ref, exc, tParity)
        integer, intent(in) :: ref(nel)
        type(Excite_0_t), intent(in) :: exc
        logical, intent(in) :: tParity

        sltcnd_excit_Excite_0_t = sltcnd_0(ref, exc, tParity)
    end function

    !>  @brief
    !>      Evaluate Matrix Element for Excite_1_t.
    !>
    !>  @param[in] ref, An array of occupied orbital indices in the reference.
    !>  @param[in] exc, An excitation of type Excite_1_t.
    !>  @param[in] tParity, The parity of the excitation.
    HElement_t(dp) function sltcnd_excit_Excite_1_t(ref, exc, tParity)
        integer, intent(in) :: ref(nel)
        type(Excite_1_t), intent(in) :: exc
        logical, intent(in) :: tParity

        sltcnd_excit_Excite_1_t = sltcnd_1(ref, exc, tParity)
    end function

    !>  @brief
    !>      Evaluate Matrix Element for Excite_1_t.
    !>
    !>  @param[in] ref, The occupied spin orbitals of the reference.
    !>  @param[in] exc, An excitation of type Excite_1_t.
    !>  @param[in] tParity, The parity of the excitation.
    HElement_t(dp) function sltcnd_excit_SpinOrbIdx_t_Excite_1_t(ref, exc, tParity)
        type(SpinOrbIdx_t), intent(in) :: ref
        type(Excite_1_t), intent(in) :: exc
        logical, intent(in) :: tParity

        sltcnd_excit_SpinOrbIdx_t_Excite_1_t = sltcnd_1(ref%idx, exc, tParity)
    end function

    !>  @brief
    !>      Evaluate Matrix Element for Excite_2_t.
    !>
    !>  @param[in] ref, An array of occupied orbital indices in the reference.
    !>  @param[in] exc, An excitation of type Excite_2_t.
    !>  @param[in] tParity, The parity of the excitation.
    HElement_t(dp) function sltcnd_excit_Excite_2_t(ref, exc, tParity)
        integer, intent(in) :: ref(nel)
        type(Excite_2_t), intent(in) :: exc
        logical, intent(in) :: tParity

        sltcnd_excit_Excite_2_t = sltcnd_2(ref, exc, tParity)
    end function

    !>  @brief
    !>      Evaluate Matrix Element for Excite_2_t.
    !>
    !>  @param[in] ref, The occupied spin orbitals of the reference.
    !>  @param[in] exc, An excitation of type Excite_2_t.
    !>  @param[in] tParity, The parity of the excitation.
    HElement_t(dp) function sltcnd_excit_SpinOrbIdx_t_Excite_2_t(ref, exc, tParity)
        type(SpinOrbIdx_t), intent(in) :: ref
        type(Excite_2_t), intent(in) :: exc
        logical, intent(in) :: tParity

        sltcnd_excit_SpinOrbIdx_t_Excite_2_t = sltcnd_2(ref%idx, exc, tParity)
    end function

    !>  @brief
    !>      Evaluate Matrix Element for Excite_3_t.
    !>
    !>  @param[in] exc, An excitation of type Excite_3_t.
    !>  @param[in] tParity, The parity of the excitation.
    HElement_t(dp) function sltcnd_excit_Excite_3_t(nI, exc, tParity)
        integer, intent(in) :: nI(nEl)
        type(Excite_3_t), intent(in) :: exc
        logical, intent(in) :: tParity

        sltcnd_excit_Excite_3_t = sltcnd_3(nI, exc, tParity)
    end function

    !>  @brief
    !>      Excitations further than max_excit_rank should return 0
    !>
    !>  @param[in] exc
    HElement_t(dp) function sltcnd_excit_Excite_Further_t(nI, exc, tParity)
        integer, intent(in) :: nI(nEl)
        type(Excite_Further_t), intent(in) :: exc
        logical, intent(in) :: tParity
        ! Only the type, not the value of this variable is used.
        @:unused_var(nI, exc, tParity)

        sltcnd_excit_Excite_Further_t = h_cast(0.0_dp)
    end function

end module
