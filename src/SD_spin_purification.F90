#include "macros.h"


module SD_spin_purification_mod
    use constants, only: n_int, dp, int64
    use growing_buffers, only: buffer_int_1D_t
    use util_mod, only: stop_all, operator(.isclose.)
    implicit none

    private
    public :: S2_expval, spin_momentum, spin_q_num, get_open_shell

contains

    ! TODO(@Oskar): Can it be complex?

    pure function S2_expval(nI, nJ) result(res)
        !! Evaluates \(< D_i | S^2 | D_j > )\
        !!
        !! Only the singly occupied (open-shell) orbitals
        !! are relevant for the evaluation.
        !! Is nonzero only, if the spin projection of bra and ket are the same.
        !! Uses the second quantisation expression for the spin operator
        !! with \(S^2 = S_z ( S_z - 1) + S_{+} S_{-})\.
        !!
        !! Since $S_z$ is basically a particle counting operator,
        !! the first summand only appears in the diagonal, i.e. `nI == nJ`.
        !!
        !! The second summand is nonzero only, if $D_i$ and $D_j$
        !! differ not at all, or if they differ by exactly one spin exchange.
        !! In the former case it evaluates to the number of $\alpha$ electrons,
        !! in the latter case it is always one.
        integer, intent(in) :: nI(:), nJ(:)
            !! The bra and ket Slater determinants in nI format.
        real(dp) :: res
            !! The matrix element.
        integer, allocatable :: oS_nI(:), oS_nJ(:)
            !! The open shell of the respective input determinants.
            !! (Can be empty!)
        logical, allocatable :: alpha_I(:), alpha_J(:)
            !! If nI and nJ have their open shell at the same spatial orbitals,
            !! this array contains true if the i-th open-shell orbital is alpha
            !! occupied.
            !! [1, 2, 3, 6] get_open_shell-> [3, 6] %2-> [1, 0] ==0-> [False, True]
            !! [1, 2, 4, 5] get_open_shell-> [4, 5] %2-> [0, 1] ==0-> [True, False]
        real(dp) :: spin_proj_term, raise_low_term
        character(*), parameter :: this_routine = 'S2_expval'

        ASSERT(size(nI) == size(nJ))

        oS_nI = get_open_shell(nI)
        oS_nJ = get_open_shell(nJ)

        res = 0.0_dp
        if (size(oS_nI) /= size(oS_nJ)) return
        if (size(oS_nI) == 0 .or. size(oS_nJ) == 0) return
        ! Are the same spatial orbitals open shell?
        ! < [a, a, 0, 0] | S^2 | [0, 0, a, a] >    ==  0
        if (any(oS_nI + mod(oS_nI, 2) /= oS_nJ + mod(oS_nJ, 2))) return

        ! This explicit allocation should not be necessary because of automatic allocation,
        ! but it does not work in ifort 18.0.x
        ! Remove the statement, if support of ifort <= 18 is dropped.
        allocate(alpha_I(size(oS_nI)), alpha_J(size(oS_nJ)))

        alpha_I = mod(oS_nI, 2) == 0
        alpha_J = mod(oS_nJ, 2) == 0

        ASSERT(count(alpha_I) == count(alpha_J))

        if (all(nI == nJ)) then
            block
                real(dp) :: s_z
                ! N_a = count(mod(nI, 2) == 0)
                ! N_b = size(nI) - N_a
                ! s_z = (N_a - N_b) / 2._dp
                ! s_z = (2 * N_a - size(nI)) / 2._dp
                s_z = (2 * count(alpha_I) - size(alpha_I)) / 2._dp
                spin_proj_term = s_z * (s_z - 1_dp)
            end block
        else
            spin_proj_term = 0._dp
        end if

        ! There is one exchange, rest is equal
        if (count(alpha_I .neqv. alpha_J) == 2) then
            raise_low_term = 1._dp
        ! Everything is equal
        else if (all(alpha_I .eqv. alpha_J)) then
            raise_low_term = real(count(alpha_I), dp)
        else
            raise_low_term = 0._dp
        end if

        res = spin_proj_term + raise_low_term
    end function

    pure function get_open_shell(nI) result(res)
        !! Return only the SOMOs.
        integer, intent(in) :: nI(:)
            !! Determinant in nI format.
        integer, allocatable :: res(:)
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
        !! Solves the equation \(X = \sqrt(s \cdot (s + 1))\)
        !! for $s$.
        real(dp), intent(in) :: spin_momentum
            !! The spin angular momentum.
        real(dp) :: res
        character(*), parameter :: this_routine = 'spin_q_num'
        ASSERT(spin_momentum >= 0._dp)
        res = -0.5_dp + sqrt(0.25_dp + spin_momentum**2)
    end function
end module SD_spin_purification_mod
