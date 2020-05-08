#include "macros.h"

! a small module with functions needed for the unit-tests
module unit_test_helpers

    use constants, only: dp, EPS, n_int, bits_n_int, maxExcit, int64, iout, lenof_sign

    use lattice_mod, only: get_helement_lattice, lattice

    use Determinants, only: get_helement, write_det

    use SystemData, only: t_lattice_model, nOccAlpha, nOccBeta, &
                          trans_corr_param_2body, omega, nel, nBasis, &
                          arr, brr, nBasis, bhub

    use fcimcdata, only: excit_gen_store_type

    use util_mod, only: binary_search, choose, operator(.div.), operator(.isclose.)

    use sltcnd_mod, only: dyn_sltcnd_excit

    use bit_reps, only: decode_bit_det, extract_sign, get_weight

    use bit_rep_data, only: niftot, nifd

    use semi_stoch_procs, only: global_most_populated_states

    use DetBitOps, only: ilut_lt, ilut_gt, EncodeBitDet, findbitexcitlevel, &
                         count_open_orbs

    use orb_idx_mod, only: SpinOrbIdx_t, new_write_det => write_det

    use excitation_types, only: Excitation_t, get_excitation

    use sort_mod, only: sort

    use matrix_util, only: matrix_exponential, blas_matmul

    use CalcData, only: PgenUnitTestSpec_t

    use ras, only: sort_orbitals

    implicit none

    private
    public :: run_excit_gen_tester, batch_run_excit_gen_tester

    abstract interface
        subroutine generate_all_excits_t(nI, n_excits, det_list)
            use SystemData, only: nel
            use constants, only: n_int
            integer, intent(in) :: nI(nel)
            integer, intent(out) :: n_excits
            integer(n_int), intent(out), allocatable :: det_list(:,:)
        end subroutine generate_all_excits_t

        real(dp) function calc_pgen_t(det_I, ilutI, exc)
            import :: dp, SpinOrbIdx_t, Excitation_t, n_int, NifTot
            type(SpinOrbIdx_t), intent(in) :: det_I
            integer(n_int), intent(in) :: ilutI(0:NIfTot)
            class(Excitation_t), intent(in) :: exc
        end function calc_pgen_t

        logical function print_predicate_t(det_I, exc, pgen_diagnostic)
            import :: dp, SpinOrbIdx_t, Excitation_t
            type(SpinOrbIdx_t), intent(in) :: det_I
            class(Excitation_t), intent(in) :: exc
            real(dp), intent(in) :: pgen_diagnostic
        end function
    end interface
contains

    subroutine setup_arr_brr(in_lat)
        class(lattice), intent(in) :: in_lat

        integer :: i

        if (associated(arr)) deallocate(arr)
        allocate(arr(nBasis,2))
        if (associated(brr)) deallocate(brr)
        allocate(brr(nBasis))

        brr = [(i, i = 1, nBasis)]
        arr = 0.0_dp

        do i = 1, nbasis
            arr(i,:) = bhub * in_lat%dispersion_rel_orb(get_spatial(i))
        end do

        call sort(arr(1:nBasis,1), brr(1:nBasis), nskip = 2)
        call sort(arr(2:nBasis,1), brr(2:nBasis), nskip = 2)
!
!         print *, "arr: "
!         do i = 1, nBasis
!             print *, arr(i,:)
!         end do
!         print *, "brr: "
!         do i = 1, nBasis
!             print *, brr(i)
!         end do

    end subroutine setup_arr_brr

    function create_hamiltonian(list_nI) result(hamil)
        ! quite specific hamiltonian creation for my tests..
        integer, intent(in) :: list_nI(:,:)
        HElement_t(dp) :: hamil(size(list_nI,2),size(list_nI,2))

        integer :: i, j

        hamil = h_cast(0.0_dp)

        do i = 1, size(list_nI,2)
            do j = 1, size(list_nI,2)
                hamil(i,j) = get_helement_lattice(list_nI(:,j),list_nI(:,i))
            end do
        end do

    end function create_hamiltonian

    function create_spin_dependent_hopping(list_nI, spin_opt) result(hamil)
        ! function to create a spin-dependent hopping matrix for the
        ! exact tests of the spin-dependent hoppint transcorrelation
        ! is no spin (1 alpha, -1 beta) is given then alpha hopping is the
        ! default
        integer, intent(in) :: list_nI(:,:)
        integer, intent(in), optional :: spin_opt
        HElement_t(dp) :: hamil(size(list_nI,2),size(list_nI,2))

        integer :: i, j, spin, ex(2,maxExcit), ic
        integer(n_int) :: ilutI(0:NifTot), ilutJ(0:niftot)
        logical :: tpar

        hamil = h_cast(0.0_dp)

        if (present(spin_opt)) then
            spin = spin_opt
        else
            spin = 1
        end if

        do i = 1, size(list_nI,2)
            call EncodeBitDet(list_nI(:,i), ilutI)
            do j = 1, size(list_nI,2)
                call EncodeBitDet(list_nI(:,j), ilutJ)

                ic = findbitexcitlevel(ilutI,ilutJ)

                if (ic /= 1) cycle

                ex(1,1) = 1
                call GetBitExcitation(ilutI,ilutJ,ex,tpar)

                if (.not. same_spin(ex(1,1),ex(2,1))) cycle

                if (get_spin_pn(ex(1,1)) == spin) then
                    hamil(i,j) = get_helement_lattice(list_nI(:,j),list_nI(:,i))
                end if

            end do
        end do

    end function create_spin_dependent_hopping

    function create_hamiltonian_old(list_nI) result(hamil)
        ! try to also create the hamiltonian with the old implementation..
        ! although i think there needs to be more setup done..
        integer, intent(in) :: list_nI(:,:)
        HElement_t(dp) :: hamil(size(list_nI,2),size(list_nI,2))

        integer :: i, j

        t_lattice_model = .false.
        do i = 1, size(list_nI,2)
            do j = 1, size(list_nI,2)
                hamil(i,j) = get_helement(list_nI(:,j),list_nI(:,i))
            end do
        end do
        t_lattice_model = .true.

    end function create_hamiltonian_old

    subroutine find_degeneracies(e_values, ind, pairs)
        ! find the indices of degenerate eigenvalues
        ! ind will have as many rows as degenerate eigenvalues exist
        ! and the columns will be the maximum number of degeneracy + 1
        ! since in the first column the number of degenerate eigenvalues are
        ! stored!
        ! it assumes the eigenvalues are sorted!!
        ! in pairs the paired indices of the degenerate eigenvalue are stored!
        real(dp), intent(in) :: e_values(:)
        integer, intent(out), allocatable :: ind(:,:), pairs(:,:)

        integer :: i,j,tmp_ind(size(e_values),size(e_values)+1), e_ind
        integer :: max_val

        tmp_ind = 0
        e_ind = 1
        i = 1

        do while (i < size(e_values) .and. e_ind < size(e_values))
            j = 0
            do while(e_ind + j <= size(e_values))
                if (abs(e_values(e_ind) - e_values(e_ind+j)) < 10.e-8) then
                    tmp_ind(i,j+2) = e_ind+j
                    j = j + 1
                else
                    exit
                end if
            end do
            tmp_ind(i,1) = j
            i = i + 1
            e_ind = e_ind + j
        end do

        ! deal with end-value specifically
        if (e_ind == size(e_values)) then
            tmp_ind(i,1) = 1
            tmp_ind(i,2) = e_ind
        end if

        max_val = maxval(tmp_ind(:,1))+1
        allocate(ind(i-1,max_val), source = tmp_ind(1:i-1,1:max_val))

        if (max_val == 2) then
            ! if no degeneracies
            allocate(pairs(size(e_values),1))
            pairs = 0
            return
        end if
        allocate(pairs(size(e_values),max_val-2))
        pairs = 0
        ! do it in a stupid way and reuse the created array ind
        do i = 1, size(ind,1)
            if (ind(i,1) > 1) then
                do j = 2, ind(i,1) + 1
                    pairs(ind(i,j),:) = pack(ind(i,2:ind(i,1)+1), &
                        ind(i,2:ind(i,1)+1) /= ind(i,j))
                end do
            end if
        end do

    end subroutine find_degeneracies

    function similarity_transform(H, t_mat_opt) result(trans_H)
        HElement_t(dp), intent(in) :: H(:,:)
        HElement_t(dp), intent(in), optional :: t_mat_opt(:,:)
        real(dp) :: trans_H(size(H,1),size(H,2))

        HElement_t(dp) :: t_mat(size(H,1),size(H,2))

        if (present(t_mat_opt)) then
            t_mat = t_mat_opt
        else
            ! otherwise assume the on-site correlation factor is used
            t_mat = get_tranformation_matrix(H, nOccAlpha*nOccBeta)
        end if

        trans_H = blas_matmul(blas_matmul(matrix_exponential(-t_mat), H), matrix_exponential(t_mat))

    end function similarity_transform



    function get_tranformation_matrix(hamil, n_pairs) result(t_matrix)
        ! n_pairs is actually also a global system dependent quantitiy..
        ! which actually might be helpful.. but input it here!
        HElement_t(dp), intent(in) :: hamil(:,:)
        integer, intent(in) :: n_pairs
        real(dp) :: t_matrix(size(hamil,1),size(hamil,2))

        integer :: i, j

        t_matrix = 0.0_dp

        do i = 1, size(hamil,1)
            do j = 1, size(hamil,1)
                if (i == j) then
                    t_matrix(i,i) = n_pairs
                else
                    if (abs(hamil(i,j)) > EPS) then
                        t_matrix(i,j) = sign(1.0_dp, real(hamil(i,j),dp))
                    end if
                end if
            end do
        end do

        t_matrix = trans_corr_param_2body/omega * t_matrix

    end function get_tranformation_matrix

    function create_all_spin_flips(nI_in) result(spin_flips)
        ! takes a given spin-configuration in nI representation and
        ! creates all possible states with flipped spin and same ms
        integer, intent(in) :: nI_in(:)
        integer, allocatable :: spin_flips(:,:)

        integer :: nI(size(nI_in)), nJ(size(nI_in))
        integer :: num, n_open, ms, n_states, i, j, k, n_found
        integer(n_int) :: ilutI(0:niftot), ilutJ(0:niftot)
        integer, allocatable :: open_shells(:)
        integer(n_int), allocatable :: ilut_list(:,:)

        nI = nI_in

        num = size(nI)

        ! it is not necessary to have nI sorted at input, so do that now
        call sort_orbitals(nI)

        call EncodeBitDet(nI, ilutI)

        n_open = count_open_orbs(ilutI)
        ms = sum(get_spin_pn(nI))

        ! the number of possible spin distributions:
        n_states = int(choose(n_open, n_open/2 + ms))

        allocate(spin_flips(num,n_states))
        spin_flips = 0
        ! the first det will be the original one
        spin_flips(:,1) = nI

        allocate(ilut_list(0:niftot,n_states))
        ilut_list = 0_n_int
        ilut_list(:,1) = ilutI

        n_found = 1
        ! i have a brute force solution for now:
        ! loop over all the possible states, starting with the original,
        ! since we need to account for multiple flips in this way
        ! i dont need to take the last one..
        do i = 1, n_states - 1
            ! here determine the positions of the open-shell orbitals
            open_shells = find_open_shell_indices(spin_flips(:,i))

            ! and now flip all possible open-shell pairs
            do j = 1, n_open - 1
                do k = j + 1, n_open
                    ! cycle if same orbital or parallel spin
!                     if (j == k) cycle
                    if (same_spin(spin_flips(open_shells(j),i), spin_flips(open_shells(k),i))) cycle

                    nj = spin_flips(:,i)
                    if (is_beta(spin_flips(open_shells(j),i))) then
                        ! then we know (j) is beta and (k) is alpha
                        nJ(open_shells(j)) = &
                            get_alpha(spin_flips(open_shells(j),i))
                        nJ(open_shells(k)) = &
                            get_beta(spin_flips(open_shells(k),i))
                    else
                        ! otherwise (j) is alpha and (k) is beta
                        nJ(open_shells(j)) = &
                            get_beta(spin_flips(open_shells(j),i))
                        nJ(open_shells(k)) = &
                            get_alpha(spin_flips(open_shells(k),i))
                    end if

                    ! if this determinant is not yet in the spin-flip
                    ! list, put it in
                    call EncodeBitDet(nJ, ilutJ)
                    if (.not. is_in_list_ilut(ilutJ, n_found, ilut_list)) then
                        n_found = n_found + 1
                        ilut_list(:,n_found) = ilutJ
                        spin_flips(:,n_found) = nJ
                    end if
                end do
            end do
        end do

!         call print_matrix(transpose(spin_flips))

    end function create_all_spin_flips

    logical function is_in_list_ilut(tgt_ilut, n_states, ilut_list_in, t_sorted_opt)
        ! interfaced function for finding a ilut rep. in a ilut list
        integer(n_int), intent(in) :: tgt_ilut(0:niftot)
        integer, intent(in) :: n_states
        integer(n_int), intent(in) :: ilut_list_in(0:niftot,n_states)
        logical, intent(in), optional :: t_sorted_opt

        logical :: t_sorted
        integer :: pos
        integer(n_int) :: ilut_list(0:niftot,n_states)

        if (present(t_sorted_opt)) then
            t_sorted = t_sorted_opt
        else
            ! default it is believed that the list is NOT sorted
            t_sorted = .false.
        end if

        ilut_list = ilut_list_in

        if (.not. t_sorted) then
            call sort(ilut_list, ilut_lt, ilut_gt)
        end if

        pos = binary_search(ilut_list, tgt_ilut, nifd+1)

        if (pos > 0) then
            is_in_list_ilut = .true.
        else
            is_in_list_ilut = .false.
        end if

    end function is_in_list_ilut

    function find_open_shell_indices(nI) result(open_shells)
        ! find the indices of the open-shell orbitals of a nI rep.
        integer, intent(in) :: nI(:)
        integer, allocatable :: open_shells(:)

        integer(n_int) :: ilutI(0:niftot)
        integer :: n_open, i, cnt

        call EncodeBitDet(nI, ilutI)

        n_open = count_open_orbs(ilutI)

        allocate(open_shells(n_open))
        open_shells = 0

        cnt = 0
        do i = 1, size(nI)
            if (.not. IsDoub(ilutI,nI(i))) then
                cnt = cnt + 1
                open_shells(cnt) = i
            end if
        end do

    end function find_open_shell_indices


!>  @brief
!>      Test if an excitation generator generates all and only expected states
!>      with the correct pgen.
!>
!>  @author Werner Dobrautz, Oskar Weser
!>
!>  @details
!>  The pgen_diagnostic is given by
!>    \f[\sum_i \frac{1}{p_i N}\f]
!>  and should be roughly one.
!>  All problematic states with respect to that diagnostic are
!>  printed in the end with their pgen_diagnostic,
!>  excitation type (ic), matrix element, and pgen.
!>  @param[in] excit_gen, An excitation generator.
!>  @param[in] excit_gen_name, The name of the excitation generator.
!>  @param[in] opt_nI, An optional reference state.
!>  @param[in] opt_n_iters, An optional number of iterations. Defaults to 100000.
!>  @param[in] gen_all_excits, An optional subroutine to generate all states
!>      that can be reached from the reference state.
!>  @param[in] calc_pgen, An optional function that calculates the pgen
!>      for a given reference an excitation. If a state is never generated,
!>      the pgen cannot be taken from the excitation generator.
!>      Adds an additional column to the output table.
!>  @param[in] print_predicate, An optional predicate function.
!>      If it returns true, an entry in the final table is printed.
!>      By default all states with a pgen_diagnostic that deviates with more than 5\,\%
!>      are printed.
    subroutine run_excit_gen_tester(excit_gen, excit_gen_name, opt_nI, opt_n_iters, &
            gen_all_excits, calc_pgen, print_predicate)
        use procedure_pointers, only: generate_excitation_t
        use util_mod, only: binary_search

        procedure(generate_excitation_t) :: excit_gen
        integer, intent(in), optional :: opt_nI(nel), opt_n_iters
        procedure(generate_all_excits_t), optional :: gen_all_excits
        procedure(calc_pgen_t), optional :: calc_pgen
        procedure(print_predicate_t), optional :: print_predicate
        character(*), intent(in) :: excit_gen_name
        character(*), parameter :: this_routine = "run_excit_gen_tester"

        integer :: i, nI(nel), n_iters
        integer :: default_n_iters = 100000
        integer :: default_n_dets = 1

        integer(n_int) :: ilut(0:niftot), tgt_ilut(0:niftot)
        integer :: nJ(nel), n_excits, ex(2,maxExcit), ic, ex_flag, i_unused = 0
        type(excit_gen_store_type) :: store
        logical :: tPar, found_all
        real(dp) :: pgen, contrib
        real(dp), allocatable :: pgen_list(:)
        HElement_t(dp) :: hel
        integer(n_int), allocatable :: det_list(:,:)
        real(dp), allocatable :: contrib_list(:)
        logical, allocatable :: generated_list(:)
        integer :: n_dets, n_generated, pos

        procedure(print_predicate_t), pointer :: print_predicate_

        ASSERT(nel > 0)
        ! and also nbasis and stuff..
        ASSERT(nbasis > 0)
        ASSERT(nel <= nbasis)

        if (present(print_predicate)) then
            print_predicate_ => print_predicate
        else
            print_predicate_ => default_predicate
        end if

        ! use some default values if not provided:
        ! nel must be set!
        if (present(opt_nI)) then
            nI = opt_nI
        else
            ! use HF-det as default
            nI = [(i, i = 1, nel)]
        end if

        if (present(opt_n_iters)) then
            n_iters = opt_n_iters
        else
            n_iters = default_n_iters
        end if

        ! i have to rewrite this routine, to a part which
        ! creates all excitations and another who runs on possibly
        ! multiple excitations!
!         if (present(opt_n_dets)) then
!             n_dets = opt_n_dets
!         else
            n_dets = default_n_dets
!         end if

        ! the problem here is now.. we want to calulate all the possible
        ! excitations.. i would need a general routine, which does that
        ! with only the hamiltonian knowledge.. this is not really there..
        ! i should have a routine:
        ! for this special setup, which is tested..

        ! for some special systems we should provide a routine to
        ! calculate all the excitations (hubbard, UEG eg.)
        if (present(gen_all_excits)) then
            call gen_all_excits(nI, n_excits, det_list)
        else
            call gen_all_excits_default(nI, n_excits, det_list)
        end if

        ! Lists to keep track of things
        allocate(generated_list(n_excits), source=.false.)
        allocate(contrib_list(n_excits), source=0.0_dp)
        allocate(pgen_list(n_excits), source=0.0_dp)

        n_generated = 0
        contrib = 0.0_dp

        n_dets = min(n_dets, n_excits)

        print *, "---------------------------------"
        print *, "testing: ", excit_gen_name
        print *, "for ", n_dets, " determinants"
        print *, " and ", n_iters, " iterations "

        call EncodeBitDet(nI, ilut)
        print *, "for starting determinant: ", nI

        write(iout, *) ! linebreak
        write(iout, '(A)') 'Progressbar'

        ! call this below now for the number of specified determinants
        ! (also use excitations of the first inputted, to be really
        !   consistent)

        block
            integer(int64) :: i, L
            L = n_iters .div. 100
            do i = 1, int(n_iters, kind=int64)
                if (mod(i, L) == 0_int64) then
                    write(iout, '(A)', advance='no') '#'
                    flush(iout)
                end if
                call excit_gen(nI, ilut, nJ, tgt_ilut, ex_flag, ic, ex, tpar, pgen, &
                            hel, store)

                if (nJ(1) == 0) cycle
                call EncodeBitDet(nJ, tgt_ilut)
                pos = binary_search(det_list, tgt_ilut, nifd+1)
                if (pos < 0) then
                    print *, "nJ: ", nJ
                    print *, "ilutJ:", tgt_ilut
                    call stop_all(this_routine, 'Unexpected determinant generated')
                else
                    generated_list(pos) = .true.
                    n_generated = n_generated + 1

                    contrib = contrib + 1.0_dp / pgen
                    contrib_list(pos) = contrib_list(pos) + 1.0_dp / pgen
                    pgen_list(pos) = pgen
                end if
            end do
            write(iout, *) ! linebreak
        end block


        print *, n_generated, " dets generated in ", n_iters, " iterations "
        print *, 100.0_dp * (n_iters - n_generated) / real(n_iters,dp), "% abortion rate"
        print *, "Averaged contribution: ", contrib / real(n_excits*n_iters,dp)

        block
            type(SpinOrbIdx_t) :: det_I
            real(dp) :: pgen_diagnostic
            class(Excitation_t), allocatable :: exc
            integer :: nJ(nEl)
            logical :: tParity


            write(iout, *) "=================================="
            write(iout, *) "Problematic contribution List: "
            write(iout, *) "=================================="
            write(iout, '("|       Determinant         |   Sum{1 / pgen} / n_iter |  ic |    <psi_I H psi_J >        |    pgen    |")', advance='no')
            if (present(calc_pgen)) write(iout, '("   calc_pgen |")', advance='no')
            write(iout, *) ! linebreak

            do i = 1, n_excits
                pgen_diagnostic = contrib_list(i) / real(n_iters, dp)
                ic = findbitexcitlevel(ilut, det_list(:, i))
                call decode_bit_det(nJ, det_list(:, i))
                call get_excitation(nI, nJ, ic, exc, tParity)

                if (print_predicate_(SpinOrbIdx_t(nI), exc, pgen_diagnostic)) then
                    call write_det (iout, nJ, .false.)
                    write(iout, '("|"F10.5"|"I2"|"F10.5"|"F15.10"|")', advance='no') &
                        pgen_diagnostic, ic, get_helement(nI, nJ),  pgen_list(i)
                    if (present(calc_pgen)) then
                        write(iout, '(F15.10"|")', advance='no') &
                            calc_pgen(SpinOrbIdx_t(nI), det_list(:, i), exc)
                    end if
                    write(iout, *)
                end if
            end do
            write(iout, *) "=================================="
            write(iout, *) ! linebreak
        end block

        contains

        logical function default_predicate(det_I, exc, pgen_diagnostic)
            type(SpinOrbIdx_t), intent(in) :: det_I
            class(Excitation_t), intent(in) :: exc
            real(dp), intent(in) :: pgen_diagnostic
            default_predicate = &
                pgen_diagnostic <= 0.95_dp .or. 1.05_dp <= pgen_diagnostic &
                .and. .not. (abs(dyn_sltcnd_excit(det_I, exc)) .isclose. 0.0_dp)
        end function

    end subroutine run_excit_gen_tester

    subroutine batch_run_excit_gen_tester(pgen_unit_test_spec)

        ! TODO(@Oskar): Remove
        use FciMCData, only: TotWalkers

        type(PgenUnitTestSpec_t), intent(in) :: pgen_unit_test_spec

        integer :: i
        type(SpinOrbIdx_t), allocatable :: largest_walkers(:)

        associate(n_most_populated => pgen_unit_test_spec%n_most_populated)
        block
            integer(n_int), allocatable :: ilut_largest_walkers(:, :)
            integer :: counter

            allocate(ilut_largest_walkers(0:NIfTot, n_most_populated), source=0_n_int)
            call global_most_populated_states(n_most_populated, ilut_largest_walkers)

            count_non_zeros: do i = 1, n_most_populated
                if (get_weight(ilut_largest_walkers(:, i)) <= 1.0e-7_dp) exit
            end do count_non_zeros
            counter = i - 1

            allocate(largest_walkers(counter))

            convert_iluts: do i = 1, counter
                ! Unfortunately there are no class methods, that can be called from the type.
                largest_walkers(i) = largest_walkers(1)%from_ilut(ilut_largest_walkers(:, i))
            end do convert_iluts
        end block
        end associate

        do i = 1, size(largest_walkers)
            call new_write_det(largest_walkers(i))
        end do
    end subroutine

    subroutine create_hilbert_space(nI, n_states, state_list_ni, state_list_ilut, &
            gen_all_excits_opt)
        ! a really basic routine, which creates the whole hilbert space based
        ! on a input determinant and other quantities, like symmetry sectors,
        ! already set outside the routine. for now this is specifically
        ! implemented for the k-space hubbard model, since i still need to
        ! test the transcorrelated approach there!
        integer, intent(in) :: nI(nel)
        integer, intent(out) :: n_states
        integer, intent(out), allocatable :: state_list_ni(:,:)
        integer(n_int), intent(out), allocatable :: state_list_ilut(:,:)
        procedure(generate_all_excits_t), optional :: gen_all_excits_opt
        character(*), parameter :: this_routine = "create_hilbert_space"

        procedure(generate_all_excits_t), pointer :: gen_all_excits
        integer(n_int), allocatable :: excit_list(:,:), temp_list_ilut(:,:)
        integer, allocatable :: temp_list_ni(:,:)
        integer :: n_excits, n_total, tmp_n_states, cnt, i, j, pos
        integer(n_int) :: ilutI(0:niftot)

        ! determine the type of system by the gen_all_excits routine

        if (present(gen_all_excits_opt)) then
            gen_all_excits => gen_all_excits_opt
        else
            gen_all_excits => gen_all_excits_default
        end if

        ! estimate the total number of excitations
        n_total = int(choose(nBasis/2, nOccAlpha) * choose(nBasis/2,nOccBeta))

        n_states = 1
        allocate(temp_list_ilut(0:niftot, n_total))
        allocate(temp_list_ni(nel, n_total))

        call EncodeBitDet(nI, ilutI)

        temp_list_ilut(:,1) = ilutI
        temp_list_ni(:,1) = nI

        tmp_n_states = 1
        ! thats a really inefficient way to do it:
        ! think of smth better at some point!
        do while (.true.)

            ! and i need a counter, which counts the number of added
            ! excitations to the whole list.. i guess if this is 0
            ! the whole hilbert space is reached..
            cnt = 0

            ! i need a temporary loop variable
            do i = 1, tmp_n_states
                call gen_all_excits(temp_list_ni(:,i), n_excits, excit_list)

                ! now i have to check if those states are already in the list
                do j = 1, n_excits

                    pos = binary_search(temp_list_ilut(:,1:(tmp_n_states + cnt)), &
                        excit_list(:,j), nifd+1)

                    ! if not yet found:
                    if (pos < 0) then
                        ! insert it at the correct place
                        ! does - pos give me the correct place then?
                        pos = -pos
                        ! lets try.. and temp_list is always big enough i think..
                        ! first move
                        temp_list_ilut(:,(pos+1):tmp_n_states+cnt+1) = &
                            temp_list_ilut(:,pos:(tmp_n_states+cnt))

                        temp_list_ni(:,(pos+1):(tmp_n_states+cnt+1)) = &
                            temp_list_ni(:,pos:(tmp_n_states+cnt))
                        ! then insert
                        temp_list_ilut(:,pos) = excit_list(:,j)

                        call decode_bit_det(temp_list_ni(:,pos), excit_list(:,j))

                        ! and increase the number of state counter
                        cnt = cnt + 1
                    else
                        ! if already found i do not need to do anything i
                        ! guess..
                    end if
                end do

            end do
            tmp_n_states = tmp_n_states + cnt

            ! and somehow i need an exit criteria, if we found all the states..
            if (cnt == 0) exit
        end do

        n_states = tmp_n_states

        ! it should be already sorted or?? i think so..
        ! or does binary_search not indicate the position
        allocate(state_list_ni(nel,n_states), source = temp_list_ni(:,1:n_states))
        allocate(state_list_ilut(0:niftot,n_states), source = temp_list_ilut(:,1:n_states))

    end subroutine create_hilbert_space

    subroutine gen_all_excits_default(nI, n_excits, det_list)
        use SymExcit3, only: CountExcitations3, GenExcitations3
        integer, intent(in) :: nI(nel)
        integer, intent(out) :: n_excits
        integer(n_int), intent(out), allocatable :: det_list(:,:)
        character(*), parameter :: this_routine = "gen_all_excits_default"

        integer :: n_singles, n_doubles, n_dets, ex(2,maxExcit), ex_flag
        integer :: nJ(nel)
        logical :: tpar, found_all
        integer(n_int) :: ilut(0:niftot)

        n_excits = -1

        call EncodeBitDet(nI, ilut)

        ! for reference in the "normal" case it looks like that:
        call CountExcitations3(nI, 2, n_singles, n_doubles)

        n_excits = n_singles + n_doubles

        print *, "n_singles: ", n_singles
        print *, "n_doubles: ", n_doubles

        allocate(det_list(0:niftot,n_excits))
        n_dets = 0
        found_all = .false.
        ex = 0
        ex_flag = 2
        call GenExcitations3 (nI, ilut, nJ, ex_flag, ex, tpar, found_all, &
                              .false.)

        do while (.not. found_all)
            n_dets = n_dets + 1
            call EncodeBitDet (nJ, det_list(:,n_dets))

            call GenExcitations3 (nI, ilut, nJ, ex_flag, ex, tpar, &
                                  found_all, .false.)
        end do

        if (n_dets /= n_excits) then
            print *, "expected number of excitations: ", n_excits
            print *, "actual calculated ones: ", n_dets
            call stop_all(this_routine,"Incorrect number of excitations found")
        end if

        ! Sort the dets, so they are easy to find by binary searching
        call sort(det_list, ilut_lt, ilut_gt)

    end subroutine gen_all_excits_default


end module unit_test_helpers
