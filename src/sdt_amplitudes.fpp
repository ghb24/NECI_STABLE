#include "macros.h"
#:include "macros.fpph"
#:include "algorithms.fpph"
#:set CI_amplitudes = ['singles_t', 'doubles_t', 'triples_t']

! Module to collect and average the CI coefficients
module sdt_amplitudes

    use bit_reps, only: extract_sign, decode_bit_det, encode_sign, niftot, nifd
    use constants, only: dp, lenof_sign, n_int, int64, stdout
    use DetBitOps, only: get_bit_excitmat, EncodeBitDet, GetBitExcitation
    use util_mod, only: near_zero, stop_all, lex_leq
    use FciMCData, only: TotWalkers, iLutRef, CurrentDets, AllNoatHF, projedet, &
                         ll_node
    use hash, only: hash_table_lookup, add_hash_table_entry, init_hash_table, &
                    clear_hash_table
    use LoggingData, only: n_store_ci_level, n_iter_ci_coeff
    use SystemData, only: nel, nbasis, symmax
    use Parallel_neci, only: iProcIndex, MPIcollection
    use MPI_wrapper, only: root
    use sort_mod, only: sort
    use fortran_strings, only: str

    better_implicit_none
    private
    public :: init_ciCoeff, store_ci_coeff, output_ci_coeff

    type, abstract :: CI_coefficients_t
    end type

    type, extends(CI_coefficients_t) :: singles_t
        real(dp) :: x
        integer :: i, a
    end type

    type, extends(CI_coefficients_t) :: doubles_t
        real(dp) :: x
        integer :: i, a, j, b
    end type

    type, extends(CI_coefficients_t) :: triples_t
        real(dp) :: x
        integer :: i, a, j, b, k, c
    end type

    integer(n_int), allocatable :: ciCoeff_storage(:, :), root_ciCoeff_storage(:, :)
    integer :: first_free_entry, nCyc, root_first_free_entry
    type(ll_node), pointer :: hash_table_ciCoeff(:)
    integer(n_int), allocatable  :: totEntCoeff(:, :)

    interface sorting
        #:for CI_exct_level in CI_amplitudes
            module procedure sorting_${CI_exct_level}$
        #:endfor
    end interface

    interface write_ci_coeff
        #:for CI_exct_level in CI_amplitudes
            module procedure write_ci_coeff_${CI_exct_level}$
        #:endfor
    end interface

contains

    subroutine init_ciCoeff
        integer, parameter :: hash_table_ciCoeff_size = 500000
        nCyc = 0
        first_free_entry = 0
        allocate (hash_table_ciCoeff(hash_table_ciCoeff_size))
        call init_hash_table(hash_table_ciCoeff)
        allocate (ciCoeff_storage(0:NIfTot, hash_table_ciCoeff_size))
        allocate (totEntCoeff(n_store_ci_level, 2))
    end subroutine init_ciCoeff


    subroutine store_ci_coeff
        integer :: ic, ex(2, 4), nIEx(nel)
        integer(int64) :: i
        real(dp) :: sign_tmp(lenof_sign)
        nCyc = nCyc + 1

        ! loop through all occupied determinants
        do i = 1, TotWalkers
            ! definition of the max ic as input for get_bit_excitmat
            ic = n_store_ci_level + 1
            ! extraction of the excitation level from every determinant
            call get_bit_excitmat(iLutRef(:, 1), CurrentDets(:, i), ex, ic)
            if (ic <= n_store_ci_level) then
                call extract_sign(CurrentDets(:, i), sign_tmp)
                call decode_bit_det(nIEx, CurrentDets(:, i))
                call cache_sign(sign_tmp, nIEx)
            end if
        end do
    end subroutine store_ci_coeff

    ! it updates the CI coeffs storage list
    subroutine cache_sign(sgn, nIEx)
        integer, intent(in) :: nIEx(nel)
        real(dp), intent(in) :: sgn(lenof_sign)

        integer :: hash_value, ind
        real(dp) :: sign_tmp(lenof_sign)
        integer(n_int) :: ilut(0:NIfTot)
        logical :: tSuccess

        ! encode the determinant into bit representation (ilut)
        call EncodeBitDet(nIEx, ilut)
        call hash_table_lookup(nIEx, ilut, NIfD, hash_table_ciCoeff, &
                               ciCoeff_storage, ind, hash_value, tSuccess)
        ! tSuccess is true when the coeff is found in the hash_table; so it gets updated
        if (tSuccess) then
            call extract_sign(ciCoeff_storage(:, ind), sign_tmp)
            sign_tmp = sign_tmp + sgn
            call encode_sign(ciCoeff_storage(:, ind), sign_tmp)

            ! tSuccess is false, then add a new entry to the CI coeffs storage
        else
            ! it counts the number of entries of different CI coeffs
            first_free_entry = first_free_entry + 1
            ! encode the determinant into bit representation (ilut)
            call EncodeBitDet(nIEx, ilut)
            ! store the encoded determinant in ciCoeff_storage
            ciCoeff_storage(:, first_free_entry) = ilut
            ! store the sign in ciCoeff_storage
            call encode_sign(ciCoeff_storage(:, first_free_entry), sgn)
            ! create a new hashtable entry
            call add_hash_table_entry(hash_table_ciCoeff, &
                                      first_free_entry, hash_value)
        end if
    end subroutine cache_sign


    ! output of the CI coefficients collection
    subroutine output_ci_coeff

        if (iProcIndex == root) then
            write (stdout, *) ''
            write (stdout, *) '=============== CI coefficients collection ==============='
            write (stdout, "(A45,I10)") 'Maximum excitation level of the CI coeff. : ', n_store_ci_level
            write (stdout, "(A45,I10)") 'Number of iterations set for average      : ', n_iter_ci_coeff
            if (nCyc /= n_iter_ci_coeff) then
                write (stdout, "(A45,I10)") ' -> actual iterations used for average    : ', nCyc
            end if
        end if

        ! it gathers all the CI coeff from different processes in one single array (root_ciCoeff_storage)
        call MPIcollection(NIfTot, first_free_entry, ciCoeff_storage, root_first_free_entry, root_ciCoeff_storage)

        if (iProcIndex == root) then
            call print_averaged_ci_coeff()

            write (stdout, *) 'Sorting CI coefficients...'
            call molpro_ci_coeff()

            write (stdout, *) '-> CI coefficients written in ASCII files ci_coeff_*'
            write (stdout, *) '=========================================================='
            write (stdout, *) ''
        end if

        call fin_ciCoeff()

    end subroutine output_ci_coeff

    ! it prints averaged CI coeffs collected during the calcualtion
    subroutine print_averaged_ci_coeff
        integer  :: i, ic, ex(2, 4), iCI, unit_CIav
        real(dp) :: sign_tmp(lenof_sign), ref_coef
        logical  :: tPar

        totEntCoeff(:, :) = 0
        iCILoop0: do iCI = 0, n_store_ci_level
            if (iCI /= 0) open (newunit=unit_CIav, file='ci_coeff_'//str(iCI)//'_av', status='replace')
            ! loop over the total entries of CI coefficients
            do i = 1, root_first_free_entry
                ic = n_store_ci_level + 1
                ! gets the excitation level of the CI coefficient
                call get_bit_excitmat(iLutRef(:, 1), root_ciCoeff_storage(:, i), ex, ic)
                if (iCI == ic) then
                    ! gets the value of the CI coefficient (i.e. the number of walkers)
                    call extract_sign(root_ciCoeff_storage(:, i), sign_tmp)
                    ex(1, 1) = ic
                    ! gets the sign of the CI coef (tPar=true if odd number of permutations)
                    call GetBitExcitation(ilutRef(:, 1), root_ciCoeff_storage(:, i), ex, tPar)
                    if (tPar) sign_tmp = -sign_tmp

                    select case (iCI) ! writing averaged CI coefficients
                    case (0)  ! reference
                        write (stdout, "(A45,F14.3)") 'Instantaneous number of walkers on HF     : ', AllNoatHF
                        ref_coef = -sign_tmp(1)
                        write (stdout, "(A45,F14.3)") 'Averaged number of walkers on HF          : ', -sign_tmp/nCyc
                        write (stdout, "(A45,I10)") 'Total entries of CI coefficients          : ', root_first_free_entry
                    case (1)  ! singles
                        totEntCoeff(iCI, 1) = totEntCoeff(iCI, 1) + 1        ! total entries for singles
                        if (.not. near_zero(sign_tmp(1))) then
                            totEntCoeff(iCI, 2) = totEntCoeff(iCI, 2) + 1    ! total entries for singles without zeros
                            write (unit_CIav, '(G20.12,2I5)') sign_tmp/ref_coef, ex(1, 1), ex(2, 1)
                        end if
                    case (2)  ! doubles
                        totEntCoeff(iCI, 1) = totEntCoeff(iCI, 1) + 1        ! total entries for doubles
                        if (.not. near_zero(sign_tmp(1))) then
                            totEntCoeff(iCI, 2) = totEntCoeff(iCI, 2) + 1    ! total entries for doubles without zeros
                            write (unit_CIav, '(G20.12,4I5)') sign_tmp/ref_coef, ex(1, 1), ex(2, 1), &
                                                                                 ex(1, 2), ex(2, 2)
                        end if
                    case (3)  ! triples
                        totEntCoeff(iCI, 1) = totEntCoeff(iCI, 1) + 1        ! total entries for triples
                        if (.not. near_zero(sign_tmp(1))) then
                            totEntCoeff(iCI, 2) = totEntCoeff(iCI, 2) + 1    ! total entries for triples without zeros
                            write (unit_CIav, '(G20.12,6I5)') sign_tmp/ref_coef, ex(1, 1), ex(2, 1), &
                                                             ex(1, 2), ex(2, 2), ex(1, 3), ex(2, 3)
                        end if
                    end select
                end if
            end do
            if (iCI /= 0) then
                close (unit_CIav)
                write (stdout, "(A31,I1,A12,I11)") 'total entries CI coeff. with ', iCI, 'excit.    :', totEntCoeff(iCI, 1)
                if (totEntCoeff(iCI, 1) /= totEntCoeff(iCI, 2)) then
                    write (stdout, "(A17,A27,I11)") '- without zeros', ':', totEntCoeff(iCI, 2)
                end if
            end if
        end do iCILoop0
    end subroutine print_averaged_ci_coeff

    ! reads the averaged CI coeff., sorts them with converted indices
    ! (OCC(alpha), OCC(beta), VIR(alpha), VIR(beta) ) and writes them
    subroutine molpro_ci_coeff
        integer  :: iCI, idxAlphaBetaOrbs(nbasis)
        class(CI_coefficients_t), allocatable :: CI_coeff(:)

        idxAlphaBetaOrbs = findAlphaBetaOrbs(symmax, nbasis)

        do iCI = 1, n_store_ci_level
            call read_ci_coeff(iCI, idxAlphaBetaOrbs, CI_coeff)
            call dyn_sort_ci_coeff(CI_coeff)
            call dyn_write_ci_coeff(iCI, CI_coeff)
        end do
    end subroutine molpro_ci_coeff

    subroutine read_ci_coeff(iCI, idxAlphaBetaOrbs, CI_coeff)
        integer, intent(in) :: iCI, idxAlphaBetaOrbs(nbasis)
        class(CI_coefficients_t), allocatable, intent(out) :: CI_coeff(:)
        integer :: h, ex(2, n_store_ci_level), signCI, unit_CIav
        integer(n_int) :: hI
        real(dp) :: x

        if (iCI == 1) then
            allocate (singles_t :: CI_coeff(totEntCoeff(iCI, 2)) )
        else if (iCI == 2) then
            allocate (doubles_t :: CI_coeff(totEntCoeff(iCI, 2)) )
        else if (iCI == 3) then
            allocate (triples_t :: CI_coeff(totEntCoeff(iCI, 2)) )
        end if

        open (newunit=unit_CIav, file='ci_coeff_'//str(iCI)//'_av', &
              status='old', action='read')
        do hI = 1_n_int, totEntCoeff(iCI, 2)
            read (unit_CIav, *) x, (ex(1, h), ex(2, h), h=1, iCI)

            call idxPreSort(iCI, idxAlphaBetaOrbs, ex, signCI)

            select type (CI_coeff)
            type is (singles_t)
                CI_coeff(hI) = singles_t(signCI*x, i=ex(1, 1), a=ex(2, 1))
            type is (doubles_t)
                CI_coeff(hI) = doubles_t(signCI*x, i=ex(1, 1), a=ex(2, 1),&
                                                   j=ex(1, 2), b=ex(2, 2))
            type is (triples_t)
                CI_coeff(hI) = triples_t(signCI*x, i=ex(1, 1), a=ex(2, 1),&
                           j=ex(1, 2), b=ex(2, 2), k=ex(1, 3), c=ex(2, 3))
            end select
        end do
        close (unit_CIav)
    end subroutine read_ci_coeff

    subroutine dyn_sort_ci_coeff(CI_coeff)
        class(CI_coefficients_t), intent(inout) :: CI_coeff(:)
        character(*), parameter :: this_routine = 'dyn_sort_ci_coeff'

        select type (CI_coeff)
        #:for CI_exct_level in CI_amplitudes
        type is(${CI_exct_level}$)
            call sorting(CI_coeff)
        #:endfor
        class default
            call stop_all(this_routine, 'Invalid CI_coefficients_t.')
        end select
    end subroutine dyn_sort_ci_coeff

    subroutine sorting_singles_t(singles)
        type(singles_t), intent(inout) :: singles(:)
        @:sort(singles_t, singles, rank=1, along=1, comp=singles_comp)
    contains
        logical elemental function singles_comp(p1, p2)
            type(singles_t), intent(in) :: p1, p2
            associate(idx_1 => [p1%i, p1%a], idx_2 => [p2%i, p2%a])
                singles_comp = lex_leq(idx_1, idx_2)
            end associate
        end function
    end subroutine sorting_singles_t

    subroutine sorting_doubles_t(doubles)
        type(doubles_t), intent(inout) :: doubles(:)
        @:sort(doubles_t, doubles, rank=1, along=1, comp=doubles_comp)
    contains
        logical elemental function doubles_comp(p1, p2)
            type(doubles_t), intent(in) :: p1, p2
            associate(idx_1 => [p1%j, p1%i, p1%b, p1%a], &
                      idx_2 => [p2%j, p2%i, p2%b, p2%a])
                doubles_comp = lex_leq(idx_1, idx_2)
            end associate
        end function
    end subroutine sorting_doubles_t

    subroutine sorting_triples_t(triples)
        type(triples_t), intent(inout) :: triples(:)
        @:sort(triples_t, triples, rank=1, along=1, comp=triples_comp)
    contains
        logical elemental function triples_comp(p1, p2)
            type(triples_t), intent(in) :: p1, p2
            associate(idx_1 => [p1%k, p1%j, p1%i, p1%a, p1%b, p1%c], &
                      idx_2 => [p2%k, p2%j, p2%i, p2%a, p2%b, p2%c])
                triples_comp = lex_leq(idx_1, idx_2)
            end associate
        end function
    end subroutine sorting_triples_t

    subroutine dyn_write_ci_coeff(iCI, CI_coeff)
        integer, intent(in) :: iCI
        class(CI_coefficients_t), allocatable, intent(inout) :: CI_coeff(:)
        integer :: unit_CIsrt

        open (newunit=unit_CIsrt, file=get_filename(CI_coeff), status='replace')
            select type (CI_coeff)
            #:for CI_exct_level in CI_amplitudes
            type is(${CI_exct_level}$)
                call write_ci_coeff(unit_CIsrt, CI_coeff)
            #:endfor
            end select
        close (unit_CIsrt)
    end subroutine dyn_write_ci_coeff

    pure function get_filename(CI_coeff) result(res)
        class(CI_coefficients_t), allocatable, intent(in) :: CI_coeff(:)
        character(:), allocatable :: res
        select type(CI_coeff)
        type is(singles_t)
            res = 'ci_coeff_1'
        type is(doubles_t)
            res = 'ci_coeff_2'
        type is(triples_t)
            res = 'ci_coeff_3'
        end select
    end function get_filename

    subroutine write_ci_coeff_singles_t(unit_CIsrt, CI_coeff)
        integer, intent(in) :: unit_CIsrt
        type(singles_t), intent(in) :: CI_coeff(:)
        integer :: h

        do h = 1, size(CI_coeff)
            write (unit_CIsrt, '(G20.12,2I5)') CI_coeff(h)%x, CI_coeff(h)%i, &
                                               CI_coeff(h)%a
        end do
    end subroutine write_ci_coeff_singles_t

    subroutine write_ci_coeff_doubles_t(unit_CIsrt, CI_coeff)
        integer, intent(in) :: unit_CIsrt
        type(doubles_t), intent(in) :: CI_coeff(:)
        integer :: h

        do h = 1, size(CI_coeff)
            write (unit_CIsrt, '(G20.12,4I5)') CI_coeff(h)%x, CI_coeff(h)%i, &
                                CI_coeff(h)%a, CI_coeff(h)%j, CI_coeff(h)%b
        end do
    end subroutine write_ci_coeff_doubles_t

    subroutine write_ci_coeff_triples_t(unit_CIsrt, CI_coeff)
        integer, intent(in) :: unit_CIsrt
        type(triples_t), intent(in) :: CI_coeff(:)
        integer :: h

        do h = 1, size(CI_coeff)
            write (unit_CIsrt, '(G20.12,6I5)') CI_coeff(h)%x, CI_coeff(h)%i, &
                                CI_coeff(h)%a, CI_coeff(h)%j, CI_coeff(h)%b, &
                                CI_coeff(h)%k, CI_coeff(h)%c
        end do
    end subroutine write_ci_coeff_triples_t

    ! it finds all the alpha/beta occ/unocc orbs in the reference
    ! determinant and saves the relative indices in an array
    function findAlphaBetaOrbs(symmax, nbasis)  result(idxAlphaBetaOrbs)
        use SymExcitDataMod, only: OrbClassCount
        use GenRandSymExcitNUMod, only: ClassCountInd
        integer, intent(in) :: symmax, nbasis
        integer :: idxAlphaBetaOrbs(nbasis)

        integer :: i, z, iSym, totEl, totOrbs
        integer :: alphaOrbs(symmax), betaOrbs(symmax), orbsSym(symmax+1)
        logical :: checkOccOrb

        orbsSym(:) = 0
        i = 1
        do iSym = 0, symmax-1
            i = i + 1
            orbsSym(i) = OrbClassCount(ClassCountInd(1, iSym, 0))*2 + orbsSym(i - 1)
        end do
        ! loop to find all the occupied alpha spin orbitals
        totEl = 0
        do iSym = 1, symmax
            alphaOrbs(iSym) = 0
            do z = 1, nel
                if (projEDet(z, 1) > orbsSym(iSym) .and. projEDet(z, 1) <= orbsSym(iSym + 1)) then
                    if (MOD(projEDet(z, 1), 2) == 1) then
                        alphaOrbs(iSym) = alphaOrbs(iSym) + 1
                        idxAlphaBetaOrbs(totEl + alphaOrbs(iSym)) = projEDet(z, 1)
                    end if
                end if
            end do
            totEl = totEl + alphaOrbs(iSym)
            ! loop to find all the occupied beta spin orbitals
            betaOrbs(iSym) = 0
            do z = 1, nel
                if (projEDet(z, 1) > orbsSym(iSym) .and. projEDet(z, 1) <= orbsSym(iSym + 1)) then
                    if (MOD(projEDet(z, 1), 2) == 0) then
                        betaOrbs(iSym) = betaOrbs(iSym) + 1
                        idxAlphaBetaOrbs(totEl + betaOrbs(iSym)) = projEDet(z, 1)
                    end if
                end if
            end do
            totEl = totEl + betaOrbs(iSym)
        end do
        if (totEl /= nel) write (stdout, *) 'WARNING: not matching number of electrons!'
        ! loop to find all the non-occupied alpha spin orbitals
        totOrbs = 0
        do iSym = 1, symmax
            alphaOrbs(iSym) = 0
            do i = orbsSym(iSym) + 1, orbsSym(iSym + 1)
                checkOccOrb = .false.
                do z = 1, nel
                    if (i == projEDet(z, 1)) checkOccOrb = .true.
                end do
                if (MOD(i, 2) == 1 .and. .not.checkOccOrb) then
                    alphaOrbs(iSym) = alphaOrbs(iSym) + 1
                    idxAlphaBetaOrbs(nel + totOrbs + alphaOrbs(iSym)) = i
                end if
            end do
            totOrbs = totOrbs + alphaOrbs(iSym)
            ! loop to find all the non-occupied beta spin orbitals
            betaOrbs(iSym) = 0
            do i = orbsSym(iSym) + 1, orbsSym(iSym + 1)
                checkOccOrb = .false.
                do z = 1, nel
                    if (i == projEDet(z, 1)) checkOccOrb = .true.
                end do
                if (MOD(i, 2) == 0 .and. .not.checkOccOrb) then
                    betaOrbs(iSym) = betaOrbs(iSym) + 1
                    idxAlphaBetaOrbs(nel + totOrbs + betaOrbs(iSym)) = i
                end if
            end do
            totOrbs = totOrbs + betaOrbs(iSym)
        end do
        if (nel + totOrbs /= nbasis) write (stdout, *) 'WARNING: not matching number of orbitals!'
    end function findAlphaBetaOrbs

    ! PreSorting routine which converts the indices into Molpro standard
    subroutine idxPreSort(iCI, idxAlphaBetaOrbs, ex, signCI)
        use util_mod, only: swap
        integer, intent(in)  :: iCI, idxAlphaBetaOrbs(nbasis)
        integer, intent(out) :: signCI
        integer, intent(inout) :: ex(2, n_store_ci_level)
        integer  :: j, k, p

        ! indices conversion
        do k = 1, iCI
            do p = 1, nel
                if (ex(1, k) == idxAlphaBetaOrbs(p)) then
                    ex(1, k) = p
                    exit
                end if
            end do
            do p = nel + 1, nbasis
                if (ex(2, k) == idxAlphaBetaOrbs(p)) then
                    ex(2, k) = p - nel
                    exit
                end if
            end do
        end do
        ! insertion sort
        signCI = 1
        do p = 1, 2
            do k = 2, iCI
                do j = k, 2, -1
                    if (ex(p, j) < ex(p, j - 1)) then
                        call swap(ex(p, j), ex(p, j - 1))
                        ! minus sign for odd number of permutations
                        signCI = -signCI
                    else
                        exit
                    end if
                end do
            end do
        end do
    end subroutine idxPreSort

    subroutine fin_ciCoeff
        call clear_hash_table(hash_table_ciCoeff)
        deallocate (hash_table_ciCoeff)
        deallocate (ciCoeff_storage)
        deallocate (totEntCoeff)
    end subroutine fin_ciCoeff

end module sdt_amplitudes
