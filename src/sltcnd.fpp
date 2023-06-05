#include "macros.h"
#:include "macros.fpph"
#:set max_excit_rank = 3
#:set excit_ranks = list(range(max_excit_rank + 1))
#:set excitations = ['Excite_{}_t'.format(i) for i in excit_ranks + ['Further']]
#:set defined_excitations = excitations[:-1]
#:set trivial_excitations = [excitations[0], excitations[-1]]
#:set non_trivial_excitations = excitations[1:-1]

#:def assert_occupation_allowed(nI, exc, assert_occupation)
#ifdef DEBUG_
    block
        use constants, only: stderr
        use util_mod, only: stop_all
        logical :: test_occupation

        if (present(${assert_occupation}$)) then
            test_occupation = ${assert_occupation}$
        else
            test_occupation = .true.
        end if

        if (test_occupation) then
            if (.not. occupation_allowed(${nI}$, ${exc}$)) then
                write(stderr, *) 'src', ${exc}$%val(1, :)
                write(stderr, *) 'tgt', ${exc}$%val(2, :)
                write(stderr, *) 'nI', ${nI}$
                call stop_all(this_routine, "Not allowed by occupation.")
            end if
        end if
    end block
#endif
#:enddef

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
                          t_mol_3_body, t_ueg_3_body, tContact, t_calc_adjoint
    use SystemData, only: nBasis!, iSpinSkip
    ! HACK - We use nBasisMax(2,3) here rather than iSpinSkip, as it appears
    !        to be more reliably set (see for example test H2O_RI)
    ! TODO: We need to sort this out so that they are consistent
    !       --> Talk to George/Alex to see what impact that might have?
    use sets_mod, only: operator(.in.), operator(.notin.), subset, disjoint, operator(.complement.)
    use constants, only: dp, n_int, maxExcit
    use UMatCache, only: GTID, UMatInd
    use IntegralsData, only: UMAT, t_use_tchint_lib
    use OneEInts, only: GetTMatEl, TMat2D
    use procedure_pointers, only: get_umat_el
    use excitation_types, only: Excitation_t, Excite_0_t, Excite_1_t, Excite_2_t, &
                                Excite_3_t, UNKNOWN, get_excitation, get_bit_excitation, &
                                create_excitation, Excite_Further_t, dyn_nI_excite, &
                                occupation_allowed, is_canonical, canonicalize
    use orb_idx_mod, only: SpinOrbIdx_t
    use DetBitOps, only: count_open_orbs, FindBitExcitLevel
    use bit_rep_data, only: NIfTot
    use LMat_mod, only: get_lmat_el, get_lmat_el_ua, external_lMat_matel
    use gen_coul_ueg_mod, only: get_contact_umat_el_3b_sp, get_contact_umat_el_3b_sap
    use SD_spin_purification_mod, only: possible_purification_methods, SD_spin_purification, &
                spin_pure_J, S2_expval_exc, nI_invariant_S2_expval_exc, ladder_op_exc
    use util_mod, only: stop_all, operator(.implies.)

    better_implicit_none
    private
    public :: initSltCndPtr, &
              nI_invariant_sltcnd_excit, sltcnd_excit, sltcnd_2_kernel, &
              dyn_sltcnd_excit_old, dyn_sltcnd_excit, &
              diagH_after_exc, sltcnd_compat, sltcnd, sltcnd_knowIC, &
              CalcFockOrbEnergy, sumfock, sltcnd_0_base, sltcnd_0_tc

    abstract interface
        HElement_t(dp) function sltcnd_0_t(nI, exc)
            import :: dp, nel, Excite_0_t
            integer, intent(in) :: nI(nel)
            type(Excite_0_t), intent(in) :: exc
        end function


        #:for rank, excite_t in zip(excit_ranks[1:], excitations[1:])
            HElement_t(dp) function sltcnd_${rank}$_t(nI, exc, tParity, assert_occupation) result(hel)
                import :: dp, nel, ${excite_t}$
                integer, intent(in) :: nI(nel)
                type(${excite_t}$), intent(in) :: exc
                logical, intent(in) :: tParity
                logical, intent(in), optional :: assert_occupation
                    !! This argument is **only** used in debug mode.
                    !! It ensures that src_i are indeed occupied and tgt_i
                    !! are unoccupied.
                    !! It is on by default.
            end function
        #:endfor

        #:for rank, excite_t in zip(excit_ranks[1 : ], defined_excitations[1 : ])
            HElement_t(dp) function diagH_after_exc_${rank}$_t(nI, E_0, exc)
                import :: dp, nEl, ${excite_t}$
                integer, intent(in) :: nI(nEl)
                HElement_t(dp), intent(in) :: E_0
                type(${excite_t}$), intent(in) :: exc
            end function
        #:endfor

        #:for rank, excite_t in zip(excit_ranks[2:], defined_excitations[2:])
            HElement_t(dp) function nI_invariant_sltcnd_${rank}$_t(exc) result(hel)
                import :: dp, nel, ${excite_t}$
                type(${excite_t}$), intent(in) :: exc
            end function
        #:endfor
    end interface

!>  @brief
!>      Evaluate Matrix Element for different excitations
!>      using the Slater-Condon rules.
!>
!>  @details
!>  This generic function uses compile time dispatch.
!>  This means that exc cannot be just of class(Excitation_t) but has
!>  to be a proper non-polymorphic type.
!>  For run time dispatch use dyn_sltcnd_excit.
!>
!>  @param[in] exc, An excitation of a subtype of Excitation_t.
    interface sltcnd_excit
        #:for rank in excit_ranks
            procedure sltcnd_${rank}$
        #:endfor
        module procedure sltcnd_excit_Excite_Further_t
        module procedure sltcnd_excit_SpinOrbIdx_t_Excite_1_t
        module procedure sltcnd_excit_SpinOrbIdx_t_Excite_2_t
    end interface


    interface diagH_after_exc
        !! Evaluate the energy of a new determinant \(D_j\) quickly.
        !!
        !! The calculation of a diagonal term of the hamiltonian,
        !! scales quadratically with the number of particles \( \mathcal{O}(N_{e}^2) \).
        !! Often we start from a determinant \(D_i\), where we know the diagonal term,
        !! and excite to a new determinant \(D_j\).
        !! Under this circumstance we can calculate the diagonal element of \(D_j\)
        !! in \( \mathcal{O}(N_{e}) \) time.
        !!
        !! In the following we will derive the necessary equations.
        !! We assume the notations and conventions of the "purple book" (Helgaker et al).
        !! The diagonal term for a determinant is given as
        !! \begin{equation*}
        !!      \langle D_i | \hat{H} | D_i \rangle
        !!   =
        !!      \sum_{I \in D_i} h_{II}
        !!       + \frac{1}{2} \sum_{I \in D_i} \sum_{J \in D_i} (g_{IIJJ} - g_{IJJI})
        !! \end{equation*}
        !!
        !! We want to calculate \( \langle D_j | \hat{H} | D_j \rangle - \langle D_i | \hat{H} | D_i \rangle \).
        !! Which we do by separately calculating the difference for the one- and two-electron term.
        !!
        !! We can rewrite the one-electron term as:
        !! \begin{equation*}
        !!          \sum_{I \in D_i} h_{II}
        !!      =
        !!          \sum_{I \in D_i \cap D_j} h_{II}
        !!              + \sum_{I \in D_i \setminus D_j} h_{II}
        !! \end{equation*}
        !! Which gives
        !! \begin{equation*}
        !!          \sum_{J \in D_j} h_{JJ}
        !!          - \sum_{I \in D_i} h_{II}
        !!      =
        !!          \sum_{J \in D_j \setminus D_i} h_{JJ}
        !!              - \sum_{I \in D_i \setminus D_j} h_{II}
        !!      =
        !!          \sum_{J \in \texttt{tgt}} h_{JJ}
        !!              - \sum_{I \in \texttt{src}} h_{II}
        !! \end{equation*}
        !!
        !!
        !! For the two electron term we define
        !! \begin{equation*}
        !!    \gamma_{IJ} = (g_{IIJJ} - g_{IJJI})
        !! \end{equation*}
        !! and note the two properties
        !! \begin{equation*}
        !!      \gamma_{IJ} = \gamma_{JI}
        !! \end{equation*}
        !! \begin{equation*}
        !!      \gamma_{II} = 0
        !! \end{equation*}
        !!
        !! We write
        !! \begin{equation*}
        !!          \frac{1}{2} \sum_{I \in D_i} \sum_{J \in D_i} \gamma_{IJ}
        !!      =
        !!      \frac{1}{2} \left[
        !!                  \sum_{I \in D_i \cap D_j} \sum_{J \in D_i \cap D_j} \gamma_{IJ}
        !!                + \sum_{I \in D_i \cap D_j} \sum_{J \in D_i \setminus D_j} \gamma_{IJ}
        !!                + \sum_{I \in D_i \setminus D_j} \sum_{J \in D_i \cap D_j} \gamma_{IJ}
        !!                + \sum_{I \in D_i \setminus D_j} \sum_{J \in D_i \setminus D_j} \gamma_{IJ}
        !!      \right]
        !! \end{equation*}
        !! \begin{equation*}
        !!      =
        !!      \frac{1}{2} \left[
        !!                  \sum_{I \in D_i \cap D_j} \sum_{J \in D_i \cap D_j} \gamma_{IJ}
        !!                + 2 \Big( \sum_{I \in D_i \cap D_j} \sum_{J \in D_i \setminus D_j} \gamma_{IJ} \Big)
        !!                + \sum_{I \in D_i \setminus D_j} \sum_{J \in D_i \setminus D_j} \gamma_{IJ}
        !!      \right]
        !! \end{equation*}
        !! In the last equality we used \( \gamma_{IJ} = \gamma_{JI} \).
        !!
        !! For the difference we get:
        !! \begin{equation*}
        !!          \frac{1}{2} \sum_{I \in D_j} \sum_{J \in D_j} \gamma_{IJ}
        !!          - \frac{1}{2} \sum_{I \in D_i} \sum_{J \in D_i} \gamma_{IJ}
        !!      =
        !!      \frac{1}{2} \left[
        !!                  2 \Big( \sum_{I \in D_i \cap D_j} \sum_{J \in D_j \setminus D_i} \gamma_{IJ} \Big)
        !!                + \sum_{I \in D_j \setminus D_i} \sum_{J \in D_j \setminus D_i} \gamma_{IJ}
        !!                - 2 \Big( \sum_{I \in D_i \cap D_j} \sum_{J \in D_i \setminus D_j} \gamma_{IJ} \Big)
        !!                - \sum_{I \in D_i \setminus D_j} \sum_{J \in D_i \setminus D_j} \gamma_{IJ}
        !!      \right]
        !! \end{equation*}
        !! \begin{equation*}
        !!      =
        !!          \sum_{I \in D_i \cap D_j} \Big(
        !!              \sum_{J \in D_j \setminus D_i} \gamma_{IJ}
        !!              - \sum_{J \in D_i \setminus D_j} \gamma_{IJ}
        !!          \Big)
        !!          + \sum_{I \in D_j \setminus D_i} \sum_{J \in D_j \setminus D_i, I < J} \gamma_{IJ}
        !!          - \sum_{I \in D_i \setminus D_j} \sum_{J \in D_i \setminus D_j, I < J} \gamma_{IJ}
        !! \end{equation*}
        !! \begin{equation*}
        !!      =
        !!          \sum_{I \in D_i \cap D_j} \Big(
        !!              \sum_{J \in \texttt{tgt}} \gamma_{IJ}
        !!              - \sum_{J \in \texttt{src}} \gamma_{IJ}
        !!          \Big)
        !!          + \sum_{I \in \texttt{tgt}} \sum_{J \in \texttt{tgt}, I < J} \gamma_{IJ}
        !!          - \sum_{I \in \texttt{src}} \sum_{J \in \texttt{src}, I < J} \gamma_{IJ}
        !! \end{equation*}
        !!
        !! In total we obtain
        !! \begin{equation*}
        !!      \langle D_j | \hat{H} | D_j \rangle - \langle D_i | \hat{H} | D_i \rangle
        !!  =
        !!        \sum_{J \in \texttt{tgt}} h_{JJ}
        !!      - \sum_{I \in \texttt{src}} h_{II}
        !!      + \sum_{I \in D_i \cap D_j} \Big(
        !!              \sum_{J \in \texttt{tgt}} \gamma_{IJ}
        !!              - \sum_{J \in \texttt{src}} \gamma_{IJ}
        !!          \Big)
        !!          + \sum_{I \in \texttt{tgt}} \sum_{J \in \texttt{tgt}, I < J} \gamma_{IJ}
        !!          - \sum_{I \in \texttt{src}} \sum_{J \in \texttt{src}, I < J} \gamma_{IJ}
        !! \end{equation*}
        #:for rank in excit_ranks[1:]
            procedure diagH_after_exc_${rank}$
        #:endfor
    end interface

    interface nI_invariant_sltcnd_excit
        #:for rank in excit_ranks[2:]
            procedure nI_invariant_sltcnd_${rank}$
        #:endfor
    end interface

    #:for rank in excit_ranks
        procedure(sltcnd_${rank}$_t), pointer :: sltcnd_${rank}$ => null()
        procedure(sltcnd_${rank}$_t), pointer :: nonadjoint_sltcnd_${rank}$ => null()
    #:endfor
    #:for rank in excit_ranks[2:]
        procedure(nI_invariant_sltcnd_${rank}$_t), pointer :: nI_invariant_sltcnd_${rank}$ => null()
    #:endfor
    #:for rank in excit_ranks[1:]
        procedure(diagH_after_exc_${rank}$_t), pointer :: diagH_after_exc_${rank}$ => null()
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
                nI_invariant_sltcnd_3 => nI_invariant_sltcnd_3_base
                sltcnd_3 => sltcnd_3_use_nI_invariant
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

                sltcnd_1 => sltcnd_1_base
                nI_invariant_sltcnd_2 => nI_invariant_sltcnd_2_purify_spin
                sltcnd_2 => sltcnd_2_use_nI_invariant
                nI_invariant_sltcnd_3 => nI_invariant_sltcnd_3_base
                sltcnd_3 => sltcnd_3_use_nI_invariant

            else
                sltcnd_0 => sltcnd_0_base
                #:for rank in excit_ranks[1:]
                    diagH_after_exc_${rank}$ => diagH_after_exc_${rank}$_base
                #:endfor
                sltcnd_1 => sltcnd_1_base
                nI_invariant_sltcnd_2 => nI_invariant_sltcnd_2_base
                sltcnd_2 => sltcnd_2_use_nI_invariant
                nI_invariant_sltcnd_3 => nI_invariant_sltcnd_3_base
                sltcnd_3 => sltcnd_3_use_nI_invariant
            end if

        end if

        if (t_calc_adjoint) then ! invert all matrix element calls
            #:for rank in excit_ranks
                nonadjoint_sltcnd_${rank}$ => sltcnd_${rank}$
                sltcnd_${rank}$ => adjoint_sltcnd_${rank}$
            #:endfor
        end if
    end subroutine initSltCndPtr

    ! We have to define this wrapper because
    ! function pointers cannot be elemental.
    ! This means that, there has to be wrapper, if we want elemental
    ! functions.
    !
    ! We have to define the wrapper here in this module,
    ! since we want to give the compiler
    ! the option to inline it.
    ! After all it will be run in the innermost loops.
    HElement_t(dp) elemental function get_2el(src1, tgt1, src2, tgt2)
        !! Return the two-electron integral.
        integer, intent(in) :: src1, tgt1, src2, tgt2
            !! The index conventions can be seen
            !! [here](https://www2.fkf.mpg.de/alavi/neci/master/page/02_dev_doc/15_index_conventions.html)
        get_2el = get_umat_el(src1, tgt1, src2, tgt2)
    end function

    #:for rank, excite_t in zip(excit_ranks[1:], defined_excitations[1:])
    HElement_t(dp) function adjoint_sltcnd_${rank}$(nI, ex, tSign, assert_occupation) result(hel)
        !! returns the adjoint sltcnd of the given rank: ${rank}$
        integer, intent(in) :: nI(nel)
        type(${excite_t}$), intent(in) :: ex
        logical, intent(in) :: tSign
        logical, intent(in), optional :: assert_occupation
        integer :: nJ(nel)
        type(${excite_t}$) :: adj_exc
        routine_name("adjoint_sltcnd_${rank}$")
        ! reverse excitation matrix and pass it to a new excitation object
        @:assert_occupation_allowed(nI, ex, assert_occupation)
        adj_exc%val(1, :) = ex%val(2, :)
        adj_exc%val(2, :) = ex%val(1, :)
        nJ = dyn_nI_excite(nI, ex)
        hel = nonadjoint_sltcnd_${rank}$(nJ, adj_exc, tSign)
#ifdef CMPLX_
        hel = conjg(hel)
#endif
    end function adjoint_sltcnd_${rank}$
    #:endfor

    HElement_t(dp) function adjoint_sltcnd_0(nI, ex) result(hel)
        !! Returns the adjoint for the diagonal element \( \langle D_i | \hat{H} | D_i \rangle \)
        integer, intent(in) :: nI(nel)
        type(Excite_0_t), intent(in) :: ex
        hel = nonadjoint_sltcnd_0(nI, ex)
#ifdef CMPLX_
        hel = conjg(hel)
#endif
    end function adjoint_sltcnd_0

!>  @brief
!>      Evaluate Matrix Element for different excitations
!>      using the Slater-Condon rules.
!>
!>  @details
!>  This generic function uses run time dispatch.
!>  This means that exc can be any subtype of class(Excitation_t).
!>  For performance reason it is advised to use sltcnd_excit,
!>  if the actual type is known at compile time.
!>
!>  @param[in] ref, The reference determinant as array of occupied orbital indices.
!>  @param[in] exc, An excitation of type excitation_t.
!>  @param[in] tParity, The parity of the excitation.
    function dyn_sltcnd_excit(ref, exc, tParity, assert_occupation) result(hel)
        integer, intent(in) :: ref(nel)
        class(Excitation_t), intent(in) :: exc
        logical, intent(in) :: tParity
        logical, intent(in), optional :: assert_occupation
        HElement_t(dp) :: hel
        character(*), parameter :: this_routine = 'dyn_sltcnd_excit'

        ! The compiler has to statically know, of what type exc is.
        select type (exc)
        #:for Excitation_t in trivial_excitations
        type is (${Excitation_t}$)
            block ! This block is just a necessary workaround for ifort18
                hel = sltcnd_excit(ref, exc)
            end block
        #:endfor
        #:for Excitation_t in non_trivial_excitations
        type is (${Excitation_t}$)
            block ! This block is just a necessary workaround for ifort18
                hel = sltcnd_excit(ref, exc, tParity, assert_occupation)
            end block
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

        if (.not. (IC /= 0 .implies. (present(ex) .and. present(tParity)))) then
            call stop_all(this_routine, "ex and tParity must be provided to &
                          &sltcnd_excit for all IC /= 0")
        end if
        hel = dyn_sltcnd_excit(nI, create_excitation(IC, ex) , tParity)
    end function

    function sltcnd_compat(nI, nJ, IC) result(hel)
        integer, intent(in) :: nI(nel), nJ(nel), IC
        HElement_t(dp) :: hel

        class(Excitation_t), allocatable :: exc
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

        class(Excitation_t), allocatable :: exc
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

    pure function CalcFockOrbEnergy(Orb, HFDet) result(hel)
        ! This calculates the orbital fock energy from
        ! the one- and two- electron integrals. This
        ! requires a knowledge of the HF determinant.
        !In: Orbital (Spin orbital notation)
        !In: HFDet (HF Determinant)
        integer, intent(in) :: HFDet(nel), Orb
        integer :: idHF(NEl), idOrb, j
        HElement_t(dp) :: hel

        ! Obtain the spatial rather than spin indices if required
        idOrb = gtID(Orb)
        idHF = gtID(HFDet)

        !GetTMATEl works with spin orbitals
        hel = GetTMATEl(Orb, Orb)
        ! Sum in the two electron contributions.
        do j = 1, nel
            hel = hel + get_umat_el(idOrb, idHF(j), idOrb, idHF(j))
        end do

        ! Exchange contribution only considered if tExch set.
        ! This is only separated from the above loop to keep "if (tExch)" out
        ! of the tight loop for efficiency.
        do j = 1, nel
            ! Exchange contribution is zero if I,J are alpha/beta
            if (tReltvy .or. (G1(Orb)%Ms == G1(HFDet(j))%Ms)) then
                hel = hel - get_umat_el(idOrb, idHF(j), idHF(j), idOrb)
            end if
        end do

    end function CalcFockOrbEnergy

    pure function SumFock(nI, HFDet) result(hel)

        ! This just calculates the sum of the Fock energies
        ! by considering the one-electron integrals and
        ! the double-counting contribution
        ! to the diagonal matrix elements. This is subtracted from
        ! the sum of the fock energies to calculate diagonal
        ! matrix elements, or added to the sum of the 1-electron
        ! integrals. The HF determinant needs to be supplied.

        integer, intent(in) :: nI(nel), HFDet(nel)
        HElement_t(dp) :: hel
        integer :: i

        hel = h_cast(0._dp)
        do i = 1, nEl
            hel = hel + CalcFockOrbEnergy(nI(i), HFDet)
        end do
    end function SumFock

    pure function sltcnd_0_base(nI, exc) result(hel)
        ! Calculate the  by the SlaterCondon Rules when the two
        ! determinants are the same (so we only need to specify one).
        integer, intent(in) :: nI(nel)
        type(Excite_0_t), intent(in) :: exc
        HElement_t(dp) :: hel, hel_sing, hel_doub, hel_tmp
        integer :: id(nel), i, j

        @:unused_var(exc)

        ! Sum in the one electron integrals (KE --> TMAT)
        hel_sing = sum(GetTMATEl(nI, nI))

        ! Obtain the spatial rather than spin indices if required
        id = gtID(nI)

        ! Sum in the two electron contributions.
        hel_doub = h_cast(0._dp)
        do i = 1, nel - 1
            hel_doub = hel_doub + sum(get_2el(id(i), id(i + 1 :), id(i), id(i + 1 :)))
        end do

        ! Exchange contribution only considered if tExch set.
        ! This is only separated from the above loop to keep "if (tExch)" out
        ! of the tight loop for efficiency.
        hel_tmp = h_cast(0._dp)
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

    function sltcnd_1_base(nI, ex, tSign, assert_occupation) result(hel)
        ! Calculate the  by the Slater-Condon Rules when the two
        ! determinants differ by one orbital exactly.
        integer, intent(in) :: nI(nel)
        type(Excite_1_t), intent(in) :: ex
        logical, intent(in) :: tSign
        logical, intent(in), optional :: assert_occupation
        HElement_t(dp) :: hel
        debug_function_name("sltcnd_1_base")

        @:assert_occupation_allowed(nI, ex, assert_occupation)

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


    ! dummy function for 3-body matrix elements without tc
    function nI_invariant_sltcnd_2_base(ex) result(hel)
        type(Excite_2_t), intent(in) :: ex
        HElement_t(dp) :: hel
        ! Only non-zero contributions if Ms preserved in each term (consider
        ! physical notation).
        hel = sltcnd_2_kernel(ex)
    end function

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

    function sltcnd_0_tc(nI, exc) result(hel)
        integer, intent(in) :: nI(nel)
        type(Excite_0_t), intent(in) :: exc
        HElement_t(dp) :: hel
        integer :: i, j, k
        integer :: dummy(1,0)

        ! get the diagonal matrix element up to 2nd order
        hel = sltcnd_0_base(nI, exc)
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

    function sltcnd_1_tc(nI, ex, tSign, assert_occupation) result(hel)
        integer, intent(in) :: nI(nel)
        type(Excite_1_t), intent(in) :: ex
        logical, intent(in) :: tSign
        logical, intent(in), optional :: assert_occupation
        HElement_t(dp) :: hel
        debug_function_name("sltcnd_1_tc")
        integer :: i, j
        @:assert_occupation_allowed(nI, ex, assert_occupation)

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

    function sltcnd_2_tc(nI, exc, tSign, assert_occupation) result(hel)
        integer, intent(in) :: nI(nel)
        type(Excite_2_t), intent(in) :: exc
        logical, intent(in) :: tSign
        logical, intent(in), optional :: assert_occupation
        HElement_t(dp) :: hel
        debug_function_name("sltcnd_2_tc")
        integer :: i
        @:assert_occupation_allowed(nI, exc, assert_occupation)

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

    function sltcnd_3_tc(nI, ex, tSign, assert_occupation) result(hel)
        integer, intent(in) :: nI(nel)
        type(Excite_3_t), intent(in) :: ex
        logical, intent(in) :: tSign
        logical, intent(in), optional :: assert_occupation
        HElement_t(dp) :: hel
        debug_function_name("sltcnd_3_tc")
        integer :: dummy(1)
        @:assert_occupation_allowed(nI, ex, assert_occupation)

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
    function nI_invariant_sltcnd_3_base(ex) result(hel)
        type(Excite_3_t), intent(in) :: ex
        HElement_t(dp) :: hel
        @:unused_var(ex)
        hel = 0
    end function

    #:for rank, excite_t in zip(excit_ranks[2:], defined_excitations[2:])
        function sltcnd_${rank}$_use_nI_invariant(nI, ex, tSign, assert_occupation) result(hel)
            integer, intent(in) :: nI(nel)
            type(${excite_t}$), intent(in) :: ex
            logical, intent(in) :: tSign
            logical, intent(in), optional :: assert_occupation
            debug_function_name("sltcnd_${rank}$_use_nI_invariant")
            HElement_t(dp) :: hel
            @:assert_occupation_allowed(nI, ex, assert_occupation)
            hel = nI_invariant_sltcnd_excit(ex)
            if (tSign) hel = -hel
        end function
    #:endfor

    !------------------------------------------------------------------------------------------!
    !      slater condon rules for ultracold atoms
    !------------------------------------------------------------------------------------------!

    pure function sltcnd_0_base_ua(nI, exc) result(hel)
        ! Calculate the  by the SlaterCondon Rules when the two
        ! determinants are the same (so we only need to specify one).
        integer, intent(in) :: nI(nel)
        type(Excite_0_t), intent(in) :: exc
        HElement_t(dp) :: hel

        HElement_t(dp) :: hel_sing, hel_doub, hel_tmp
        integer :: id(nel), i, j, idN, idX

        @:unused_var(exc)

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

    function sltcnd_1_base_ua(nI, ex, tSign, assert_occupation) result(hel)
        ! Calculate the  by the Slater-Condon Rules when the two
        ! determinants differ by one orbital exactly.
        integer, intent(in) :: nI(nel)
        type(Excite_1_t), intent(in) :: ex
        logical, intent(in) :: tSign
        logical, intent(in), optional :: assert_occupation
        HElement_t(dp) :: hel
        debug_function_name("sltcnd_1_base_ua")
        @:assert_occupation_allowed(nI, ex, assert_occupation)

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

    function sltcnd_2_base_ua(nI, ex, tSign, assert_occupation) result(hel)
        ! Calculate the  by the Slater-Condon Rules when the two
        ! determinants differ by two orbitals exactly (the simplest case).
        integer, intent(in) :: nI(nel)
        type(Excite_2_t), intent(in) :: ex
        logical, intent(in) :: tSign
        logical, intent(in), optional :: assert_occupation
        HElement_t(dp) :: hel
        debug_function_name("sltcnd_2_base_ua")

        @:assert_occupation_allowed(nI, ex, assert_occupation)

        ! Only non-zero contributions if Ms preserved in each term (consider
        ! physical notation).
        hel = sltcnd_2_kernel_ua(ex)

        if (tSign) hel = -hel
    end function sltcnd_2_base_ua

    function sltcnd_2_kernel_ua(ex) result(hel)
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

    pure function sltcnd_0_tc_ua(nI, exc) result(hel)
        integer, intent(in) :: nI(nel)
        type(Excite_0_t), intent(in) :: exc
        HElement_t(dp) :: hel
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

    function sltcnd_1_tc_ua(nI, exc, tSign, assert_occupation) result(hel)
        integer, intent(in) :: nI(nel)
        type(Excite_1_t), intent(in) :: exc
        logical, intent(in) :: tSign
        logical, intent(in), optional :: assert_occupation
        HElement_t(dp) :: hel
        routine_name("sltcnd_1_tc_ua")
        integer :: i, j
        @:assert_occupation_allowed(nI, exc, assert_occupation)

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

    function sltcnd_2_tc_ua(nI, ex, tSign, assert_occupation) result(hel)
        integer, intent(in) :: nI(nel)
        type(Excite_2_t), intent(in) :: ex
        logical, intent(in) :: tSign
        logical, intent(in), optional :: assert_occupation
        debug_function_name("sltcnd_2_tc_ua")
        HElement_t(dp) :: hel, heltc
        @:assert_occupation_allowed(nI, ex, assert_occupation)

        ! get the matrix element up to 2-body terms
        hel = sltcnd_2_kernel_ua(ex)
        ! and the 3-body term
        heltc = sltcnd_2_kernel_ua_3b(nI, ex)
        hel = hel + heltc

        ! take fermi sign into account
        if (tSign) hel = -hel
    end function sltcnd_2_tc_ua

    function sltcnd_3_tc_ua(nI, ex, tSign, assert_occupation) result(hel)
        integer, intent(in) :: nI(nEl)
        type(Excite_3_t), intent(in) :: ex
        logical, intent(in) :: tSign
        logical, intent(in), optional :: assert_occupation
        debug_function_name("sltcnd_3_tc_ua")
        HElement_t(dp) :: hel
        @:assert_occupation_allowed(nI, ex, assert_occupation)

        ! this is directly the fully symmetrized entry of the L-matrix
        associate(ex => ex%val)
            hel = get_lmat_el_ua(ex(1, 1), ex(1, 2), ex(1, 3), &
                                 ex(2, 1), ex(2, 2), ex(2, 3))
        end associate
        ! take fermi sign into account
        if (tSign) hel = -hel
    end function sltcnd_3_tc_ua

    function sltcnd_0_purify_spin_only_ladder(nI, exc) result(hel)
        integer, intent(in) :: nI(nel)
        type(Excite_0_t), intent(in) :: exc
        HElement_t(dp) :: hel
        hel = sltcnd_0_base(nI, exc) + spin_pure_J * ladder_op_exc(nI, exc)
    end function

    function sltcnd_0_purify_spin_full_s2(nI, exc) result(hel)
        integer, intent(in) :: nI(nel)
        type(Excite_0_t), intent(in) :: exc
        HElement_t(dp) :: hel
        hel = sltcnd_0_base(nI, exc) + spin_pure_J * S2_expval_exc(nI, exc)
    end function


    function nI_invariant_sltcnd_2_purify_spin(exc) result(hel)
        type(Excite_2_t), intent(in) :: exc
        HElement_t(dp) :: hel
        hel = nI_invariant_sltcnd_2_base(exc) + spin_pure_J * nI_invariant_S2_expval_exc(exc)
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
        routine_name("sltcnd_excit_SpinOrbIdx_t_Excite_1_t")

        @:ASSERT(subset(exc%val(1, :), ref%idx) .and. disjoint(exc%val(2, :), ref%idx))
        sltcnd_excit_SpinOrbIdx_t_Excite_1_t = sltcnd_1(ref%idx, exc, tParity)
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
        routine_name("sltcnd_excit_SpinOrbIdx_t_Excite_2_t")

        @:ASSERT(subset(exc%val(1, :), ref%idx) .and. disjoint(exc%val(2, :), ref%idx))
        sltcnd_excit_SpinOrbIdx_t_Excite_2_t = sltcnd_2(ref%idx, exc, tParity)
    end function

    !>  @brief
    !>      Excitations further than max_excit_rank should return 0
    !>
    !>  @param[in] exc
    HElement_t(dp) function sltcnd_excit_Excite_Further_t(nI, exc)
        integer, intent(in) :: nI(nEl)
        type(Excite_Further_t), intent(in) :: exc
        ! Only the type, not the value of this variable is used.
        @:unused_var(nI, exc)

        sltcnd_excit_Excite_Further_t = h_cast(0.0_dp)
    end function



    #:for rank, excite_t in zip(excit_ranks[1 : ], defined_excitations[1 : ])
        pure function diagH_after_exc_${rank}$_base(nI, E_0, exc) result(hel)
            integer, intent(in) :: nI(nEl)
            HElement_t(dp), intent(in) :: E_0
            type(${excite_t}$), intent(in) :: exc
            HElement_t(dp) :: hel
            routine_name("diagH_after_exc_${rank}$_base")

            HElement_t(dp) :: Delta_1, Delta_2
            integer :: i, j
            integer :: i_src

            @:pure_ASSERT(is_canonical(exc))

            associate(src => exc%val(1, :), tgt => exc%val(2, :))
            associate(id_nI => gtID(nI), id_src => gtID(src), id_tgt => gtID(tgt))
                Delta_1 = sum(GetTmatEl(tgt, tgt)) - sum(GetTMatEl(src, src))

                Delta_2 = h_cast(0._dp)

                i_src = 1

                do i = 1, size(nI)
                    if (i_src <= size(src)) then
                        if (nI(i) == src(i_src)) then
                            i_src = i_src + 1
                            cycle
                        end if
                    end if

                    do j = 1, size(tgt)
                        Delta_2 = Delta_2 &
                                    + get_2el(id_nI(i), id_tgt(j), id_nI(i), id_tgt(j)) &
                                    - get_2el(id_nI(i), id_src(j), id_nI(i), id_src(j))

                        if (G1(nI(i))%Ms == G1(tgt(j))%Ms) then
                            Delta_2 = Delta_2 &
                                        - get_2el(id_nI(i), id_tgt(j), id_tgt(j), id_nI(i))
                        end if
                        if (G1(nI(i))%Ms == G1(src(j))%Ms) then
                            Delta_2 = Delta_2 &
                                        + get_2el(id_nI(i), id_src(j), id_src(j), id_nI(i))
                        end if
                    end do
                end do

                do i = 1, size(tgt)
                    do j = i + 1, size(tgt)
                        Delta_2 = Delta_2 &
                                    + get_2el(id_tgt(i), id_tgt(j), id_tgt(i), id_tgt(j)) &
                                    - get_2el(id_src(i), id_src(j), id_src(i), id_src(j))

                        if (G1(tgt(i))%Ms == G1(tgt(j))%Ms) then
                            Delta_2 = Delta_2 &
                                        - get_2el(id_tgt(i), id_tgt(j), id_tgt(j), id_tgt(i))
                        end if
                        if (G1(src(i))%Ms == G1(src(j))%Ms) then
                            Delta_2 = Delta_2 &
                                        + get_2el(id_src(i), id_src(j), id_src(j), id_src(i))
                        end if
                    end do
                end do
            end associate
            end associate
            hel = E_0 + Delta_1 + Delta_2
        end function
    #:endfor

end module
