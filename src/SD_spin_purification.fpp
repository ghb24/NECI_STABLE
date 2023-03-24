#include "macros.h"
#:include "macros.fpph"
#:set excitations = ['Excite_0_t', 'Excite_1_t', 'Excite_2_t', 'Excite_3_t']


module SD_spin_purification_mod
    use constants, only: n_int, dp, int64
    use growing_buffers, only: buffer_int_1D_t
    use util_mod, only: stop_all, operator(.isclose.), swap, operator(.div.), &
        EnumBase_t
    use sets_mod, only: subset, disjoint
    use excitation_types, only: excitation_t, Excite_0_t, Excite_1_t, Excite_2_t, &
                                Excite_3_t, UNKNOWN, get_excitation, get_bit_excitation, &
                                create_excitation, is_canonical, occupation_allowed, canonicalize
    use orb_idx_mod, only: calc_spin_raw, get_spat
    implicit none

    type, extends(EnumBase_t) :: SD_SpinPurificationMethods_t
    end type

    type :: Possible_SD_SpinPurificationMethods_t
        type(SD_SpinPurificationMethods_t) :: &
            FULL_S2 = SD_SpinPurificationMethods_t(1), &
            ONLY_LADDER = SD_SpinPurificationMethods_t(2), &
            TRUNCATED_LADDER = SD_SpinPurificationMethods_t(3)
    end type

    type(Possible_SD_SpinPurificationMethods_t), parameter :: &
        possible_purification_methods = Possible_SD_SpinPurificationMethods_t()

    type(SD_SpinPurificationMethods_t), allocatable :: SD_spin_purification

    real(dp), allocatable :: spin_pure_J

    private
    public :: S2_expval, spin_momentum, spin_q_num, get_open_shell, &
        S2_expval_exc, dyn_S2_expval_exc, &
        nI_invariant_S2_expval_exc, &
        ladder_op_exc, dyn_ladder_op_exc, &
        nI_invariant_ladder_op_exc, &
        possible_purification_methods, SD_spin_purification, &
        spin_pure_J


    public :: old_ladder_op_exc_Excite_2_t



    interface S2_expval_exc
    #:for T in excitations
        module procedure S2_expval_exc_${T}$
    #:endfor
    end interface

    interface nI_invariant_S2_expval_exc
    #:for T in excitations[1:]
        module procedure nI_invariant_S2_expval_exc_${T}$
    #:endfor
    end interface

    interface ladder_op_exc
    #:for T in excitations
        module procedure ladder_op_exc_${T}$
    #:endfor
    end interface

    interface nI_invariant_ladder_op_exc
    #:for T in excitations[1:]
        module procedure nI_invariant_ladder_op_exc_${T}$
    #:endfor
    end interface

contains

    pure function S2_expval(nI, nJ) result(res)
        !! Evaluates \(< D_i | S^2 | D_j > \)
        !!
        !! Only the singly occupied (open-shell) orbitals
        !! are relevant for the evaluation.
        !! Is nonzero only, if the spin projection of bra and ket are the same.
        !! Uses the second quantisation expression for the spin operator
        !! with \( S^2 = S_z ( S_z - 1) + S_{+} S_{-} \).
        !!
        !! Since \(S_z\) is basically a particle counting operator,
        !! the first summand only appears in the diagonal, i.e. `all(nI == nJ)`.
        !!
        !! The second summand is nonzero only, if \( D_i \) and \( D_j \)
        !! differ not at all, or if they differ by exactly one spin exchange.
        !! In the former case it evaluates to the number of open shell \( \alpha \) electrons,
        !! in the latter case it is always one.
        integer, intent(in) :: nI(:)
            !! The bra Slater determinant in nI format.
        integer, intent(in) :: nJ(:)
            !! The ket Slater determinant in nI format.
        real(dp) :: res
            !! The matrix element.
            !! It is real even for complex `NECI`.
        integer, allocatable :: oS_nI(:), oS_nJ(:)
            !! The open shell of the respective input determinants.
            !! (Can be empty!)
        logical, allocatable :: alpha_I(:), alpha_J(:)
            !! If nI and nJ have their open shell at the same spatial orbitals,
            !! this array contains true if the i-th open-shell orbital is alpha
            !! occupied.
            !! `[1, 2, 3, 6] get_open_shell-> [3, 6] %2-> [1, 0] ==0-> [False, True]`
            !! `[1, 2, 4, 5] get_open_shell-> [4, 5] %2-> [0, 1] ==0-> [True, False]`
        character(*), parameter :: this_routine = 'S2_expval'

        res = 0.0_dp
        ASSERT(size(nI) == size(nJ))

        oS_nI = get_open_shell(nI)
        if (size(oS_nI) == 0) return
        oS_nJ = get_open_shell(nJ)
        if (size(oS_nJ) == 0) return

        if (size(oS_nI) /= size(oS_nJ)) return
        ! Are the same spatial orbitals open shell?
        ! < [a, a, 0, 0] | S^2 | [0, 0, a, a] >    ==  0
        if (any(oS_nI + mod(oS_nI, 2) /= oS_nJ + mod(oS_nJ, 2))) return

        ! This explicit allocation should not be necessary because of automatic allocation,
        ! but it does not work in ifort 18.0.x
        ! Remove the statement, if support of ifort <= 18 is dropped.
        allocate(alpha_I(size(oS_nI)), alpha_J(size(oS_nJ)))

        alpha_I = mod(oS_nI, 2) == 0
        alpha_J = mod(oS_nJ, 2) == 0

        ! We assume that inside NECI all determinants have the
        ! same spin projection.
        ASSERT(count(alpha_I) == count(alpha_J))

        ! There is one exchange, rest is equal
        if (count(alpha_I .neqv. alpha_J) == 2) then
            res = 1._dp
        ! Everything is equal
        else if (all(alpha_I .eqv. alpha_J)) then
        block
            real(dp) :: s_z
            ! N_a = count(mod(nI, 2) == 0)
            ! N_b = size(nI) - N_a
            ! s_z = (N_a - N_b) / 2._dp
            ! s_z = (2 * N_a - size(nI)) / 2._dp
            s_z = (2 * count(alpha_I) - size(alpha_I)) / 2._dp
            res = s_z * (s_z - 1_dp) + real(count(alpha_I), dp)
        end block
        else
            res = 0.0_dp
        end if
    end function

    pure function dyn_S2_expval_exc(nI, exc) result(res)
        !! Evaluates \(< D_i | S^2 | D_j > \)
        !!
        !! \( D_j \) is connected to \( D_i \) via the excitation `exc`.
        !!
        !! Note that the sign in front of the off-diagonal elements
        !! is always +1 if we assume the NECI order of spin orbitals.
        integer, intent(in) :: nI(:)
            !! The bra Slater determinant in nI format.
        class(Excitation_t), intent(in) :: exc
            !! An excitation.
        real(dp) :: res

        select type(exc)
        type is (Excite_0_t)
            res = S2_expval_exc(nI, exc)
        type is (Excite_1_t)
            res = S2_expval_exc(nI, exc)
        type is (Excite_2_t)
            res = S2_expval_exc(nI, exc)
        type is (Excite_3_t)
            res = S2_expval_exc(nI, exc)
        end select
    end function

    pure function dyn_ladder_op_exc(nI, exc) result(res)
        !! Evaluates \(< D_i | S_-S_+ | D_j > \)
        !!
        !! \( D_j \) is connected to \( D_i \) via the excitation `exc`.
        !!
        !! Note that the sign in front of the off-diagonal elements
        !! is always +1 if we assume the NECI order of spin orbitals.
        integer, intent(in) :: nI(:)
            !! The bra Slater determinant in nI format.
        class(Excitation_t), intent(in) :: exc
            !! An excitation.
        real(dp) :: res

        select type(exc)
        type is (Excite_0_t)
            res = ladder_op_exc(nI, exc)
        type is (Excite_1_t)
            res = ladder_op_exc(nI, exc)
        type is (Excite_2_t)
            res = ladder_op_exc(nI, exc)
        type is (Excite_3_t)
            res = ladder_op_exc(nI, exc)
        end select
    end function

    pure function nI_invariant_ladder_op_exc_Excite_1_t(exc) result(res)
        !! Evaluates \(< D_i | S_+S- | a^\dagger_A a_I D_i > = 0 \)
        type(Excite_1_t), intent(in) :: exc
        real(dp) :: res
            !! The matrix element is always exactly zero
        @:unused_var(exc)
        res = 0.0_dp
    end function

    pure function nI_invariant_ladder_op_exc_Excite_2_t(exc) result(res)
        !! Evaluates \(< D_i | S_+ S_- | a^\dagger_A a^\dagger_B a_I a_J D_i > = 0 \)
        type(Excite_2_t), intent(in) :: exc
        real(dp) :: res
        debug_function_name("nI_invariant_ladder_op_exc_Excite_2_t")
            !! The matrix element.
            !! It is real even for complex `NECI`.
        @:pure_ASSERT(is_canonical(exc))

        ! Only exchange excitations are non-zero.
        associate(srcs => exc%val(1, :), tgts => exc%val(2, :))
            if (calc_spin_raw(srcs(1)) /= calc_spin_raw(srcs(2)) &
                 .and. all(get_spat(srcs) == get_spat(tgts))) then
                 res = 1._dp
            else
                res = 0._dp
            end if
        end associate
    end function

    pure function nI_invariant_ladder_op_exc_Excite_3_t(exc) result(res)
        !! Evaluates \(< D_i | S_+S_- | a^\dagger_A a^\dagger_B a^\dagger_C a_I a_J a_K D_i > = 0 \)
        type(Excite_3_t), intent(in) :: exc
        real(dp) :: res
            !! The matrix element is always exactly zero
        @:unused_var(exc)
        res = 0.0_dp
    end function

    pure function ladder_op_exc_Excite_0_t(nI, exc) result(res)
        !! Evaluates \(< D_i | S_+S_- | D_i > \)
        integer, intent(in) :: nI(:)
            !! The bra Slater determinant in nI format.
        type(Excite_0_t), intent(in) :: exc
        real(dp) :: res
            !! The matrix element.
            !! It is real even for complex `NECI`.
        integer, allocatable :: oS_nI(:)
            !! The open-shell spin-orbitals of nI.
            !! Can be empty (allocated, but size == 0).
        @:unused_var(exc)
        oS_nI = get_open_shell(nI)
        res = real(count(mod(oS_nI, 2) == 0 ), dp)
    end function

    pure function ladder_op_exc_Excite_1_t(nI, exc) result(res)
        !! Evaluates \(< D_i | S_+S- | a^\dagger_A a_I D_i > = 0 \)
        integer, intent(in) :: nI(:)
            !! The bra Slater determinant in nI format.
        type(Excite_1_t), intent(in) :: exc
        real(dp) :: res
            !! The matrix element is always exactly zero
        @:unused_var(nI, exc)
        res = nI_invariant_ladder_op_exc(exc)
    end function

    pure function ladder_op_exc_Excite_2_t(nI, exc) result(res)
        !! Evaluates \(< D_i | S_+ S_- | a^\dagger_A a^\dagger_B a_I a_J D_i > = 0 \)
        integer, intent(in) :: nI(:)
            !! The bra Slater determinant in nI format.
        type(Excite_2_t), intent(in) :: exc
        real(dp) :: res
            !! The matrix element.
            !! It is real even for complex `NECI`.
        debug_function_name("ladder_op_exc_Excite_2_t")
        @:pure_ASSERT(is_canonical(exc))
        if (occupation_allowed(nI, exc)) then
            res = nI_invariant_ladder_op_exc(exc)
        else
            res = 0._dp
        end if
    end function

    pure function ladder_op_exc_Excite_3_t(nI, exc) result(res)
        !! Evaluates \(< D_i | S_+S_- | a^\dagger_A a^\dagger_B a^\dagger_C a_I a_J a_K D_i > = 0 \)
        integer, intent(in) :: nI(:)
            !! The bra Slater determinant in nI format.
        type(Excite_3_t), intent(in) :: exc
        real(dp) :: res
            !! The matrix element is always exactly zero
        @:unused_var(nI, exc)
        res = nI_invariant_ladder_op_exc(exc)
    end function

    pure function nI_invariant_S2_expval_exc_Excite_1_t(exc) result(res)
        !! Evaluates \(< D_i | S^2 | a^\dagger_A a_I D_i > = 0 \)
        type(Excite_1_t), intent(in) :: exc
        real(dp) :: res
            !! The matrix element is always exactly zero
        @:unused_var(exc)
        res = 0.0_dp
    end function

    pure function nI_invariant_S2_expval_exc_Excite_2_t(exc) result(res)
        !! Evaluates \(< D_i | S^2 | a^\dagger_A a^\dagger_B a_I a_J D_i > = 0 \)
        type(Excite_2_t), intent(in) :: exc
        real(dp) :: res
            !! The matrix element.
            !! It is real even for complex `NECI`.
        res = nI_invariant_ladder_op_exc(exc)
    end function

    pure function nI_invariant_S2_expval_exc_Excite_3_t(exc) result(res)
        !! Evaluates \(< D_i | S^2 | a^\dagger_A a^\dagger_B a^\dagger_C a_I a_J a_K D_i > = 0 \)
        type(Excite_3_t), intent(in) :: exc
        real(dp) :: res
            !! The matrix element is always exactly zero
        @:unused_var(exc)
        res = 0.0_dp
    end function

    pure function S2_expval_exc_Excite_0_t(nI, exc) result(res)
        !! Evaluates \(< D_i | S^2 | D_i > \)
        integer, intent(in) :: nI(:)
            !! The bra Slater determinant in nI format.
        type(Excite_0_t), intent(in) :: exc
        real(dp) :: res
            !! The matrix element.
            !! It is real even for complex `NECI`.
        real(dp) :: s_z
            !! The spin sprojection.
        @:unused_var(exc)
        s_z = (2 * count(mod(nI, 2) == 0) - size(nI)) / 2._dp
        res = s_z * (s_z - 1_dp) + ladder_op_exc(nI, exc)
    end function

    pure function S2_expval_exc_Excite_1_t(nI, exc) result(res)
        !! Evaluates \(< D_i | S^2 | a^\dagger_A a_I D_i > = 0 \)
        integer, intent(in) :: nI(:)
            !! The bra Slater determinant in nI format.
        type(Excite_1_t), intent(in) :: exc
        real(dp) :: res
            !! The matrix element is always exactly zero
        @:unused_var(nI, exc)
        res = 0._dp
    end function

    pure function S2_expval_exc_Excite_2_t(nI, exc) result(res)
        !! Evaluates \(< D_i | S^2 | a^\dagger_A a^\dagger_B a_I a_J D_i > = 0 \)
        integer, intent(in) :: nI(:)
            !! The bra Slater determinant in nI format.
        type(Excite_2_t), intent(in) :: exc
        real(dp) :: res
            !! The matrix element.
            !! It is real even for complex `NECI`.
        if (occupation_allowed(nI, exc)) then
            res = nI_invariant_S2_expval_exc(exc)
        else
            res = 0._dp
        end if
    end function

    pure function S2_expval_exc_Excite_3_t(nI, exc) result(res)
        !! Evaluates \(< D_i | S^2 | a^\dagger_A a^\dagger_B a^\dagger_C a_I a_J a_K D_i > = 0 \)
        integer, intent(in) :: nI(:)
            !! The bra Slater determinant in nI format.
        type(Excite_3_t), intent(in) :: exc
        real(dp) :: res
            !! The matrix element is always exactly zero
        @:unused_var(nI, exc)
        res = 0._dp
    end function

    pure function get_open_shell(nI) result(res)
        !! Return only the SOMOs.
        integer, intent(in) :: nI(:)
            !! Determinant in nI format.
        integer, allocatable :: res(:)
            !! The singly occupied orbitals.
            !! Result can be empty, i.e. allocated but `size == 0`.
        integer :: i
        type(buffer_int_1D_t) :: buffer

        ! There can be at most as many open-shell orbitals,
        ! as there are occupied orbitals.
        call buffer%init(start_size=size(nI, kind=int64))

        i = 1
        do while (i <= size(nI))
            if (i == size(nI)) then
                call buffer%push_back(nI(i))
                i = i + 1
            else if (mod(nI(i), 2) == 0) then
                call buffer%push_back(nI(i))
                i = i + 1
            else if (nI(i + 1) /= nI(i) + 1) then
                call buffer%push_back(nI(i))
                i = i + 1
            else
                i = i + 2
            end if
        end do
        call buffer%dump_reset(res)
    end function

    elemental function spin_momentum(s) result(res)
        !! Return the angular momentum for a spin quantum number s.
        real(dp), intent(in) :: s
            !! The spin quantum number. \(s = \frac{n}{2}, n \in \mathbb{N}\)
        real(dp) :: res
        character(*), parameter :: this_routine = 'spin_momentum'
        ASSERT((2._dp * s) .isclose. (real(nint(2.0_dp * s, int64), dp)))
        res = sqrt(s * (s + 1_dp))
    end function

    elemental function spin_q_num(spin_momentum) result(res)
        !! Return the spin quantum number for a given angular momentum.
        !!
        !! Solves the equation \(X = \sqrt{s \cdot (s + 1)}\)
        !! for \( s \).
        real(dp), intent(in) :: spin_momentum
            !! The spin angular momentum.
        real(dp) :: res
        character(*), parameter :: this_routine = 'spin_q_num'
        ASSERT(spin_momentum >= 0._dp)
        res = -0.5_dp + sqrt(0.25_dp + spin_momentum**2)
    end function


    pure function old_ladder_op_exc_Excite_2_t(nI, exc) result(res)
        !! Evaluates \(< D_i | S_+ S_- | a^\dagger_A a^\dagger_B a_I a_J D_i > = 0 \)
        integer, intent(in) :: nI(:)
            !! The bra Slater determinant in nI format.
        type(Excite_2_t), intent(in) :: exc
        real(dp) :: res
            !! The matrix element.
            !! It is real even for complex `NECI`.
        integer, allocatable :: oS_nI(:)
        integer :: src(2), src_spat(2), tgt_spat(2)
        logical :: negative
        oS_nI = get_open_shell(nI)
        res = 0.0_dp

        if (.not. occupation_allowed(nI, canonicalize(exc))) return

        negative = .false.
        if (size(oS_nI) /= 0) then
            src(:) = exc%val(1, :)
            ! Test if double excitation is of exchange type.
            if (count(mod(src, 2) == 0) == 1) then
                if (src(1) > src(2)) then
                    call swap(src(1), src(2))
                    negative = .true. .neqv. negative
                end if
                if (subset(src, os_nI)) then
                    tgt_spat = (exc%val(2, :) + 1) .div. 2
                    if (tgt_spat(1) > tgt_spat(2)) then
                        call swap(tgt_spat(1), tgt_spat(2))
                        negative = .true. .neqv. negative
                    end if
                    src_spat = (src + 1) .div. 2
                    if (all(src_spat == tgt_spat)) then
                        res = 1.0_dp
                    end if
                end if
            end if
        end if
        if (negative) res = -res
    end function
end module SD_spin_purification_mod
