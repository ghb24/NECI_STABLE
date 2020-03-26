#include "macros.h"

! This module contains a type and routines for defining and creating a sparse Hamiltonian.
! The type created to store this information is sparse_matrix_real. For an N-by-N matrix,
! one creates a 1d array of type sparse_matrix_real. Each element of this array stores the
! information on one single row of the matrix - the positions of the non-zero elements
! and their values, in the same order. This completley defines the matrix.

! Note that in some parts of NECI, the sparse nature is used in a different form (i.e. in
! the Lanczos code).

module sparse_arrays

    use bit_rep_data, only: NIfTot, NIfDBO, NIfD
    use bit_reps, only: decode_bit_det, nifguga
    use CalcData, only: tReadPops, t_guga_mat_eles
    use constants
    use DetBitOps, only: DetBitEq, CountBits, TestClosedShellDet
    use Determinants, only: get_helement
    use FciMCData, only: determ_space_size, determ_sizes, determ_displs, &
                         SpawnedParts, Hii, core_ham_diag
    use hphf_integrals, only: hphf_diag_helement, hphf_off_diag_helement, &
                              hphf_off_diag_helement_opt
    use MemoryManager, only: TagIntType, LogMemAlloc, LogMemDealloc
    use Parallel_neci, only: iProcIndex, nProcessors, MPIBarrier, MPIAllGatherV
    use SystemData, only: tHPHF, nel
    use global_det_data, only: set_det_diagH
    use shared_rhash, only: shared_rhash_t
    use SystemData, only: tGUGA
    use guga_excitations, only: calc_off_diag_guga_gen, actHamiltonian, &
                                calc_guga_matrix_element
    use guga_bitRepOps, only: convert_ilut_toGUGA, extract_h_element, &
                              init_csf_information
    use util_mod, only: binary_search
    use guga_data, only: tag_excitations, ExcitationInformation_t
    use guga_matrixElements, only: calcDiagMatEleGuga_nI

    implicit none
    type sparse_matrix_real
        HElement_t(dp), allocatable, dimension(:) :: elements
        integer, allocatable, dimension(:) :: positions
        integer :: num_elements
    end type sparse_matrix_real

    type sparse_matrix_int
        integer, allocatable, dimension(:) :: elements
        integer, allocatable, dimension(:) :: positions
        integer :: num_elements
    end type sparse_matrix_int

    type trial_hashtable
        ! All the states with this hash value.
        integer(n_int), allocatable, dimension(:,:) :: states
        ! The number of clashes for ths hash value.
        integer :: nclash
    end type trial_hashtable

    type core_hashtable
        ! The indices of states with this hash table.
        integer, allocatable, dimension(:) :: ind
        ! The number of clashes for this hash value.
        integer :: nclash
    end type core_hashtable

    type(sparse_matrix_real), allocatable, dimension(:) :: sparse_ham
    integer(TagIntType), allocatable, dimension(:,:) :: SparseHamilTags

    ! For quick access it is often useful to have just the diagonal elements. Note,
    ! however, that they *are* stored in sparse_ham too.
    real(dp), allocatable, dimension(:) :: hamil_diag
    integer(TagIntType) :: HDiagTag

    ! The core Hamiltonian for semi-stochastiic simulations.
    type(sparse_matrix_real), allocatable, dimension(:) :: sparse_core_ham
    integer(TagIntType), allocatable, dimension(:,:) :: SparseCoreHamilTags

    ! Stores the parities for all connected pairs of states in the core space.
    type(sparse_matrix_int), allocatable, dimension(:) :: core_connections

    type(trial_hashtable), allocatable, dimension(:) :: trial_ht
    type(trial_hashtable), allocatable, dimension(:) :: con_ht
    type(shared_rhash_t) :: core_ht

    ! --- For when using the determ-proj-approx-hamil option -----
    type(shared_rhash_t) :: var_ht
    type(sparse_matrix_real), allocatable, dimension(:) :: approx_ham

contains

    subroutine calculate_sparse_hamiltonian_non_hermitian(num_states, ilut_list)
        ! same routine as below, but for non-hermitian Hamiltonians in the
        ! case of transcorrelation
        integer, intent(in) :: num_states
        integer(n_int), intent(in) :: ilut_list(0:NIfTot, num_states)
        integer :: i, j, counter, ierr
        integer :: nI(nel), nJ(nel)
        HElement_t(dp), allocatable, dimension(:) :: hamiltonian_row
        integer, allocatable, dimension(:) :: sparse_diag_positions, sparse_row_sizes, indices
        integer(TagIntType) :: HRTag, SRTag, SDTag, ITag
        character(len=*), parameter :: t_r = "calculate_sparse_hamiltonian_non_hermitian"

        allocate(sparse_ham(num_states))
        allocate(SparseHamilTags(2, num_states))
        allocate(hamiltonian_row(num_states), stat=ierr)
        call LogMemAlloc('hamiltonian_row', num_states, 8, t_r, HRTag, ierr)
        allocate(hamil_diag(num_states), stat=ierr)
        call LogMemAlloc('hamil_diag', num_states, 8, t_r, HDiagTag, ierr)
        allocate(sparse_row_sizes(num_states), stat=ierr)
        call LogMemAlloc('sparse_row_sizes', num_states, bytes_int, t_r, SRTag, ierr)

        ! Set each element to one to count the diagonal elements straight away.
        sparse_row_sizes = 1

        if (tGUGA) then
            call stop_all(t_r, "modify get_helement for GUGA!")
        end if

        do i = 1, num_states

            hamiltonian_row = 0.0_dp

            call decode_bit_det(nI, ilut_list(:,i))

            ! we have to loop over everything in case on non-hermiticity
            do j = 1, num_states

                call decode_bit_det(nJ, ilut_list(:,j))
                if (i == j) then
                    if (tHPHF) then
                        hamiltonian_row(i) = hphf_diag_helement(nI, ilut_list(:,i))
                    else
                        hamiltonian_row(i) = get_helement(nI, nI, 0)
                    end if
                    hamil_diag(i) = hamiltonian_row(i)
                else
                    if (tHPHF) then
                        !TODO: do i need <I|H|J> or <J|H|I>?
                        hamiltonian_row(j) = hphf_off_diag_helement(&
                            nI, nJ, ilut_list(:,i), ilut_list(:,j))
                    else
                        hamiltonian_row(j) = get_helement(&
                            nI, nJ, ilut_list(:,i), ilut_list(:,j))
                    end if
                    if (abs(hamiltonian_row(j)) > EPS) then
                        ! i think in the non-hermitian i only need to update
                        ! one of the counters..
                        sparse_row_sizes(i) = sparse_row_sizes(i) + 1
                    end if
                end if
            end do

            call allocate_sparse_ham_row(sparse_ham, i, sparse_row_sizes(i), &
                "sparse_ham", SparseHamilTags(:,i))

            sparse_ham(i)%elements = 0.0_dp
            sparse_ham(i)%positions = 0
            sparse_ham(i)%num_elements = sparse_row_sizes(i)

            ! now fill in the elements, all of them
            counter = 1
            do j = 1, num_states
                if (abs(hamiltonian_row(j)) > EPS) then
                    sparse_ham(i)%positions(counter) = j
                    sparse_ham(i)%elements(counter) = hamiltonian_row(j)
                    counter = counter + 1
                end if
            end do

        end do

        deallocate(hamiltonian_row, stat=ierr)
        call LogMemDealloc(t_r, HRTag, ierr)
        deallocate(sparse_row_sizes, stat=ierr)
        call LogMemDealloc(t_r, SRTag, ierr)

    end subroutine calculate_sparse_hamiltonian_non_hermitian

    subroutine calculate_sparse_hamiltonian(num_states, ilut_list)

        integer, intent(in) :: num_states
        integer(n_int), intent(in) :: ilut_list(0:NIfTot, num_states)
        integer :: i, j, counter, ierr
        integer :: nI(nel), nJ(nel)
        HElement_t(dp), allocatable, dimension(:) :: hamiltonian_row
        integer, allocatable, dimension(:) :: sparse_diag_positions, sparse_row_sizes, indices
        integer(TagIntType) :: HRTag, SRTag, SDTag, ITag
        character(len=*), parameter :: t_r = "calculate_sparse_hamiltonian"

        integer :: pos, nexcits
        integer(n_int) :: ilutG(0:nifguga)
        integer(n_int), pointer :: excitations(:,:)
        type(ExcitationInformation_t) :: excitInfo

        allocate(sparse_ham(num_states))
        allocate(SparseHamilTags(2, num_states))
        allocate(hamiltonian_row(num_states), stat=ierr)
        call LogMemAlloc('hamiltonian_row', num_states, 8, t_r, HRTag, ierr)
        allocate(hamil_diag(num_states), stat=ierr)
        call LogMemAlloc('hamil_diag', num_states, 8, t_r, HDiagTag, ierr)
        allocate(sparse_row_sizes(num_states), stat=ierr)
        call LogMemAlloc('sparse_row_sizes', num_states, bytes_int, t_r, SRTag, ierr)
        allocate(sparse_diag_positions(num_states), stat=ierr)
        call LogMemAlloc('sparse_diag_positions', num_states, bytes_int, t_r, SDTag, ierr)
        allocate(indices(num_states), stat=ierr)
        call LogMemAlloc('indices', num_states, bytes_int, t_r, ITag, ierr)

        ! Set each element to one to count the diagonal elements straight away.
        sparse_row_sizes = 1

        ! also need to change this routine here for the guga implementation

        do i = 1, num_states

            hamiltonian_row = 0.0_dp

            call decode_bit_det(nI, ilut_list(:, i))

            ! sparse_diag_positions(i) stores the number of non-zero elements in row i of the
            ! Hamiltonian, up to and including the diagonal element.
            ! sparse_row_sizes(i) stores this number currently as all non-zero elements before
            ! the diagonal have been counted (as the Hamiltonian is symmetric).
            sparse_diag_positions(i) = sparse_row_sizes(i)

            if (tGUGA .and. (.not. t_guga_mat_eles)) then

                call convert_ilut_toGUGA(ilut_list(:,i), ilutG)

                call actHamiltonian(ilutG, excitations, nexcits)

                do j = 1, num_states

                    if (i == j) then
                        ! diag case
                        hamiltonian_row(j) = get_helement(nI, nI, 0)

                        hamil_diag(j) = hamiltonian_row(j)

                    else
                        ! off-diagonal case

                        pos = binary_search(excitations(0:nifd,1:nexcits), ilut_list(0:nifd,j))

                        if (pos > 0) then

                            hamiltonian_row(j) = extract_h_element(excitations(:,pos))

                            sparse_row_sizes(i) = sparse_row_sizes(i) + 1
                            sparse_row_sizes(j) = sparse_row_sizes(j) + 1

                        end if
                    end if
                end do
                deallocate(excitations)
                ! am i sure if i want to do that all the time???
                call LogMemDealloc(t_r, tag_excitations)

            else
                do j = i, num_states

                    call decode_bit_det(nJ, ilut_list(:, j))

                    ! If on the diagonal of the Hamiltonian.
                    if (i == j) then
                        if (tHPHF) then
                            hamiltonian_row(j) = hphf_diag_helement(nI, ilut_list(:, i))
                        else if (tGUGA) then
                            hamiltonian_row(j) = calcDiagMatEleGuga_nI(nI)
                        else
                            hamiltonian_row(j) = get_helement(nI, nJ, 0)
                        end if
                        hamil_diag(j) = hamiltonian_row(j)
                    else
                        if (tHPHF) then
                            hamiltonian_row(j) = hphf_off_diag_helement(nI, nJ, ilut_list(:, i), &
                                                                               ilut_list(:, j))
                        else if (tGUGA) then
                            call calc_guga_matrix_element(ilut_list(:,j), ilut_list(:,i), &
                                    excitInfo, hamiltonian_row(j), .true., 2)
                        else
                            hamiltonian_row(j) = get_helement(nI, nJ, ilut_list(:, i), &
                                                                     ilut_list(:, j))
                        end if
                        if (abs(hamiltonian_row(j)) > 0.0_dp) then
                            ! If element is nonzero, update the following sizes.
                            sparse_row_sizes(i) = sparse_row_sizes(i) + 1
                            sparse_row_sizes(j) = sparse_row_sizes(j) + 1
                        end if
                    end if
                end do
            end if
            ! Now we know the number of non-zero elements in this row of the Hamiltonian, so allocate it.
            call allocate_sparse_ham_row(sparse_ham, i, sparse_row_sizes(i), "sparse_ham", SparseHamilTags(:,i))

            sparse_ham(i)%elements = 0.0_dp
            sparse_ham(i)%positions = 0
            sparse_ham(i)%num_elements = sparse_row_sizes(i)

            ! Now fill in all matrix elements beyond and including the diagonal, as these are
            ! stored in hamiltonian_row.
            counter = sparse_diag_positions(i)
            do j = i, num_states
                if (abs(hamiltonian_row(j)) > 0.0_dp) then
                    sparse_ham(i)%positions(counter) = j
                    sparse_ham(i)%elements(counter) = hamiltonian_row(j)
                    counter = counter + 1
                end if
            end do
        end do

        ! At this point, sparse_ham has been allocated with the correct sizes, but only the
        ! matrix elements above and including the diagonal have been filled in. Now we must
        ! fill in the other elements. To do this, cycle through every element above the
        ! diagonal and fill in every corresponding element below it:
        indices = 1
        do i = 1, num_states
            do j = sparse_diag_positions(i) + 1, sparse_row_sizes(i)

                sparse_ham(sparse_ham(i)%positions(j))%&
                    &positions(indices(sparse_ham(i)%positions(j))) = i
                sparse_ham(sparse_ham(i)%positions(j))%&
                    &elements(indices(sparse_ham(i)%positions(j))) = h_conjg(sparse_ham(i)%elements(j))

                indices(sparse_ham(i)%positions(j)) = indices(sparse_ham(i)%positions(j)) + 1

            end do
        end do

        deallocate(hamiltonian_row, stat=ierr)
        call LogMemDealloc(t_r, HRTag, ierr)
        deallocate(sparse_row_sizes, stat=ierr)
        call LogMemDealloc(t_r, SRTag, ierr)
        deallocate(sparse_diag_positions, stat=ierr)
        call LogMemDealloc(t_r, SDTag, ierr)
        deallocate(indices, stat=ierr)
        call LogMemDealloc(t_r, ITag, ierr)

    end subroutine calculate_sparse_hamiltonian

    subroutine calculate_sparse_ham_par(num_states, ilut_list, tPrintInfo)

        integer(MPIArg), intent(in) :: num_states(0:nProcessors-1)
        integer(n_int), intent(in) :: ilut_list(0:NIfTot, num_states(iProcIndex))
        logical, intent(in) :: tPrintInfo
        integer(MPIArg) :: disps(0:nProcessors-1)
        integer :: i, j, row_size, counter, num_states_tot, ierr, bytes_required
        integer :: nI(nel), nJ(nel)
        integer(n_int), allocatable, dimension(:,:) :: temp_store
        integer(TagIntType) :: TempStoreTag, HRTag, SDTag
        HElement_t(dp), allocatable, dimension(:) :: hamiltonian_row
        character(len=*), parameter :: t_r = "calculate_sparse_ham_par"

        integer :: pos, nexcits
        integer(n_int) :: ilutG(0:nifguga)
        integer(n_int), pointer :: excitations(:,:)
        type(ExcitationInformation_t) :: excitInfo

        num_states_tot = int(sum(num_states), sizeof_int)
        disps(0) = 0
        do i = 1, nProcessors-1
            disps(i) = disps(i-1) + num_states(i-1)
        end do

        safe_realloc_e(sparse_ham, (num_states(iProcIndex)), ierr)
        safe_realloc_e(SparseHamilTags, (2, num_states(iProcIndex)), ierr)
        safe_realloc_e(hamiltonian_row, (num_states_tot), ierr)
        call LogMemAlloc('hamiltonian_row', num_states_tot, 8, t_r, HRTag, ierr)
        safe_realloc_e(hamil_diag, (num_states(iProcIndex)), ierr)
        call LogMemAlloc('hamil_diag', int(num_states(iProcIndex),sizeof_int), 8, t_r, HDiagTag, ierr)
        safe_realloc_e(temp_store, (0:NIfTot, num_states_tot), ierr)
        call LogMemAlloc('temp_store', num_states_tot*(NIfTot+1), 8, t_r, TempStoreTag, ierr)

        ! Stick together the determinants from all processors, on all processors.
        call MPIAllGatherV(ilut_list(:,1:num_states(iProcIndex)), temp_store, num_states, disps)

        ! Loop over all determinants on this processor.
        do i = 1, num_states(iProcIndex)

            call decode_bit_det(nI, ilut_list(:,i))

            row_size = 0
            hamiltonian_row = 0.0_dp

            if (tGUGA .and. (.not. t_guga_mat_eles)) then
                call convert_ilut_toGUGA(ilut_list(:,i), ilutG)

                call actHamiltonian(ilutG, excitations, nexcits)

                do j = 1, num_states_tot

                    if (DetBitEq(ilut_list(:,i), temp_store(:,j), nifdbo)) then

                        hamiltonian_row(j) = get_helement(nI, nI, 0)

                        hamil_diag(i) = hamiltonian_row(j)

                        row_size = row_size + 1

                    else

                        pos = binary_search(excitations(0:nifd,1:nexcits), temp_store(0:nifd,j))

                        if (pos > 0) then

                            hamiltonian_row(j) = extract_h_element(excitations(:,pos))

                            row_size = row_size + 1

                        end if
                    end if
                end do
                deallocate(excitations)
                call LogMemDealloc(t_r, tag_excitations)
            else
            ! Loop over all determinants on all processors.
            do j = 1, num_states_tot

                call decode_bit_det(nJ, temp_store(:,j))

                ! If on the diagonal of the Hamiltonian.
                if (DetBitEq(ilut_list(:,i), temp_store(:,j), NIfDBO)) then
                    if (tHPHF) then
                        hamiltonian_row(j) = hphf_diag_helement(nI, ilut_list(:,i))
                    else if (tGUGA) then
                        hamiltonian_row(j) = calcDiagMatEleGuga_nI(nI)
                    else
                        hamiltonian_row(j) = get_helement(nI, nJ, 0)
                    end if
                    hamil_diag(i) = hamiltonian_row(j)
                    ! Always include the diagonal elements.
                    row_size = row_size + 1
                else
                    if (tHPHF) then
                        hamiltonian_row(j) = hphf_off_diag_helement(nI, nJ, ilut_list(:,i), temp_store(:,j))
                    else if (tGUGA) then
                        call calc_guga_matrix_element(temp_store(:,j), ilut_list(:,i), &
                            excitInfo, hamiltonian_row(j), .true., 2)
                    else
                        hamiltonian_row(j) = get_helement(nI, nJ, ilut_list(:,i), temp_store(:,j))
                    end if
                    if (abs(hamiltonian_row(j)) > 0.0_dp) row_size = row_size + 1
                end if
            end do

            if (tPrintInfo) then
                if (i == 1) then
                    bytes_required = row_size*(8+bytes_int)
                    write(6,'(1x,a43)') "About to allocate first row of Hamiltonian."
                    write(6,'(1x,a40,1x,i8)') "The memory (bytes) required for this is:", bytes_required
                    write(6,'(1x,a71,1x,i7)') "The total number of determinants (and hence rows) on this processor is:", &
                                               num_states(iProcIndex)
                    write(6,'(1x,a58,1x,i7)') "The total number of determinants across all processors is:", num_states_tot
                    write(6,'(1x,a77,1x,i7)') "It is therefore expected that the total memory (MB) required will be roughly:", &
                                               num_states_tot*bytes_required/1000000
                else if (mod(i,1000) == 0) then
                    write(6,'(1x,a23,1x,i7)') "Finished computing row:", i
                end if
            end if

            ! Now we know the number of non-zero elements in this row of the Hamiltonian, so allocate it.
            call allocate_sparse_ham_row(sparse_ham, i, row_size, "sparse_ham", SparseHamilTags(:,i))

            sparse_ham(i)%elements = 0.0_dp
            sparse_ham(i)%positions = 0
            sparse_ham(i)%num_elements = row_size

            counter = 1
            do j = 1, num_states_tot
                ! If non-zero or a diagonal element.
                if (abs(hamiltonian_row(j)) > 0.0_dp .or. (j == i + disps(iProcIndex)) ) then
                    sparse_ham(i)%positions(counter) = j
                    sparse_ham(i)%elements(counter) = hamiltonian_row(j)
                    counter = counter + 1
                end if
                if (counter == row_size + 1) exit
            end do

        end if
        end do

        call MPIBarrier(ierr)

        deallocate(temp_store, stat=ierr)
        call LogMemDealloc(t_r, TempStoreTag, ierr)
        deallocate(hamiltonian_row, stat=ierr)
        call LogMemDealloc(t_r, HRTag, ierr)

    end subroutine calculate_sparse_ham_par

    subroutine calc_determ_hamil_sparse()

        use SystemData, only: t_3_body_excits,t_mol_3_body,t_ueg_transcorr
        integer :: i, j, row_size, counter, ierr
        integer :: nI(nel), nJ(nel)
        integer(n_int), allocatable, dimension(:,:) :: temp_store
        integer, allocatable :: temp_store_nI(:,:)
        integer(TagIntType) :: HRTag, TempStoreTag
        HElement_t(dp), allocatable, dimension(:) :: hamiltonian_row

        integer :: pos, nExcit
        integer(n_int) :: ilutG(0:nifguga)
        integer(n_int), pointer :: excitations(:,:)
        type(ExcitationInformation_t) :: excitInfo

        character(len=*), parameter :: t_r = "calc_determ_hamil_sparse"

        integer(n_int) :: tmp(0:NIfD)
        integer :: IC

        allocate(sparse_core_ham(determ_sizes(iProcIndex)), stat=ierr)
        allocate(SparseCoreHamilTags(2, determ_sizes(iProcIndex)))
        allocate(hamiltonian_row(determ_space_size), stat=ierr)
        call LogMemAlloc('hamiltonian_row', int(determ_space_size,sizeof_int), 8, t_r, HRTag, ierr)
        allocate(core_ham_diag(determ_sizes(iProcIndex)), stat=ierr)
        allocate(temp_store(0:NIfTot, determ_space_size), stat=ierr)
        call LogMemAlloc('temp_store', determ_space_size*(NIfTot+1), 8, t_r, TempStoreTag, ierr)
        safe_realloc_e(temp_store_nI, (nel, determ_space_size), ierr)

        ! Stick together the deterministic states from all processors, on
        ! all processors.
        ! n.b. Explicitly use 0:NIfTot, as NIfTot may not equal NIfBCast
        call MPIAllGatherV(SpawnedParts(0:NIfTot, 1:determ_sizes(iProcIndex)),&
                           temp_store, determ_sizes, determ_displs)

        do i = 1, determ_space_size
            call decode_bit_det(temp_store_nI(:,i), temp_store(:,i))
        end do

        ! Loop over all deterministic states on this processor.
        do i = 1, determ_sizes(iProcIndex)

            !call decode_bit_det(nI, SpawnedParts(:, i))
            nI = temp_store_nI(:, i + determ_displs(iProcIndex))

            row_size = 0
            hamiltonian_row = 0.0_dp

            ! here i should do something different for the guga case..
            ! and just apply the hamiltonian once, and then check if the
            ! other deterministic states are connected to nI
            ! make this optional, dependent on how i want to calculate the
            ! off-diagoanal elements for the guga case

            if (tGUGA .and. (.not. t_guga_mat_eles)) then
                call convert_ilut_toGUGA(SpawnedParts(:,i), ilutG)

                call actHamiltonian(ilutG, excitations, nExcit)

                ! then loop over j
                do j = 1, determ_space_size

                    if (all(SpawnedParts(0:nifdbo,i) == temp_store(0:nifdbo,j))) then

                        ! thats the diagonal case
                        hamiltonian_row(j) = get_helement(nI, nI, 0) - Hii

                        core_ham_diag(i) = hamiltonian_row(j)

                        if (.not. tReadPops) then
                            call set_det_diagH(i, real(hamiltonian_row(j), dp))
                        end if

                        row_size = row_size + 1

                    else
                        ! here i need the off-diagonal element

                        pos = binary_search(excitations(0:nifd,1:nExcit), temp_store(0:nifd,j))

                        if (pos > 0) then
                            hamiltonian_row(j) = extract_h_element(excitations(:,pos))
                            row_size = row_size + 1
                        end if
                    end if
                end do

                deallocate(excitations)
                call LogMemDealloc(t_r, tag_excitations)

            else
                ! Loop over all deterministic states.
                do j = 1, determ_space_size

                    !call decode_bit_det(nJ, temp_store(:,j))
                    nJ = temp_store_nI(:,j)

                    ! If on the diagonal of the Hamiltonian.
                    if (all( SpawnedParts(0:NIfDBO, i) == temp_store(0:NIfDBO, j) )) then
                        if (tHPHF) then
                            hamiltonian_row(j) = hphf_diag_helement(nI, SpawnedParts(:,i)) - Hii
                        else
                            ! for guga: the diagonal is fine, since i overwrite
                            ! that within get_helement
                            hamiltonian_row(j) = get_helement(nI, nJ, 0) - Hii
                        end if
                        core_ham_diag(i) = hamiltonian_row(j)
                        ! We calculate and store the diagonal matrix element at
                        ! this point for later access.
                        if (.not. tReadPops) &
                            call set_det_diagH(i, Real(hamiltonian_row(j), dp))
                        ! Always include the diagonal elements.
                        row_size = row_size + 1
                    else
                        if (tHPHF) then
                            hamiltonian_row(j) = hphf_off_diag_helement(nI, nJ, SpawnedParts(:,i), temp_store(:,j))
                        else if (tGUGA) then
                            ! for the off-diagonal elements i have to call the GUGA
                            ! specific function
                            ! but this is a waste.. i do not have to do that for
                            ! every nJ i could just check the list generated
                            ! by H|nI>..
                            call calc_guga_matrix_element(temp_store(:,j), SpawnedParts(:,i), &
                                excitInfo, hamiltonian_row(j), .true., 2)
                        else
                            tmp = ieor(SpawnedParts(0:NIfD,i), temp_store(0:NIfD,j))
                            tmp = iand(SpawnedParts(0:NIfD,i), tmp)
                            IC = CountBits(tmp, NIfD)

                            if (IC <= maxExcit) then
                                hamiltonian_row(j) = get_helement(nI, nJ, IC, SpawnedParts(:, i), temp_store(:, j))
                            end if

                        end if
                        if (abs(hamiltonian_row(j)) > 0.0_dp) row_size = row_size + 1
                    end if

                end do
            end if
            ! Now we know the number of non-zero elements in this row of the Hamiltonian, so allocate it.
            call allocate_sparse_ham_row(sparse_core_ham, i, row_size, "sparse_core_ham", SparseCoreHamilTags(:,i))

            sparse_core_ham(i)%elements = 0.0_dp
            sparse_core_ham(i)%positions = 0
            sparse_core_ham(i)%num_elements = row_size

            counter = 1
            do j = 1, determ_space_size
                ! If non-zero or a diagonal element.
                if (abs(hamiltonian_row(j)) > 0.0_dp .or. (j == i + determ_displs(iProcIndex)) ) then
                    sparse_core_ham(i)%positions(counter) = j
                    sparse_core_ham(i)%elements(counter) = hamiltonian_row(j)
                    counter = counter + 1
                end if
                if (counter == row_size + 1) exit
            end do

        end do

        ! Don't time the mpi_barrier call, because if we did then we wouldn't
        ! be able separate out some of the core Hamiltonian creation time from
        ! the MPIBarrier calls in the main loop.
        call MPIBarrier(ierr, tTimeIn=.false.)

        deallocate(temp_store, stat=ierr)
        call LogMemDealloc(t_r, TempStoreTag, ierr)
        deallocate(hamiltonian_row, stat=ierr)
        call LogMemDealloc(t_r, HRTag, ierr)
        deallocate(temp_store_nI, stat=ierr)

    end subroutine calc_determ_hamil_sparse

    subroutine calc_determ_hamil_sparse_hphf()

        integer :: i, j, row_size, counter, ierr
        integer :: nI(nel), nJ(nel)
        integer(n_int), allocatable, dimension(:,:) :: temp_store
        integer, allocatable :: temp_store_nI(:,:)
        integer(TagIntType) :: HRTag, TempStoreTag
        HElement_t(dp), allocatable, dimension(:) :: hamiltonian_row
        character(len=*), parameter :: t_r = "calc_determ_hamil_sparse_hphf"

        integer(n_int) :: tmp(0:NIfD)
        integer :: IC
        logical :: CS_I
        logical, allocatable :: cs(:)

        allocate(sparse_core_ham(determ_sizes(iProcIndex)), stat=ierr)
        allocate(SparseCoreHamilTags(2, determ_sizes(iProcIndex)))
        allocate(hamiltonian_row(determ_space_size), stat=ierr)
        call LogMemAlloc('hamiltonian_row', int(determ_space_size,sizeof_int), 8, t_r, HRTag, ierr)
        allocate(core_ham_diag(determ_sizes(iProcIndex)), stat=ierr)
        allocate(temp_store(0:NIfTot, determ_space_size), stat=ierr)
        call LogMemAlloc('temp_store', determ_space_size*(NIfTot+1), 8, t_r, TempStoreTag, ierr)
        allocate(temp_store_nI(nel, determ_space_size), stat=ierr)
        allocate(cs(determ_space_size), stat=ierr)

        ! Stick together the deterministic states from all processors, on
        ! all processors.
        ! n.b. Explicitly use 0:NIfTot, as NIfTot may not equal NIfBCast
        call MPIAllGatherV(SpawnedParts(0:NIfTot, 1:determ_sizes(iProcIndex)),&
                           temp_store, determ_sizes, determ_displs)

        do i = 1, determ_space_size
            call decode_bit_det(temp_store_nI(:,i), temp_store(:,i))
            cs(i) = TestClosedShellDet(temp_store(:,i))
        end do

        ! Loop over all deterministic states on this processor.
        do i = 1, determ_sizes(iProcIndex)

            !call decode_bit_det(nI, SpawnedParts(:, i))
            nI = temp_store_nI(:, i + determ_displs(iProcIndex))

            row_size = 0
            hamiltonian_row = 0.0_dp

            CS_I = cs(i + determ_displs(iProcIndex))

            ! Loop over all deterministic states.
            do j = 1, determ_space_size

                !call decode_bit_det(nJ, temp_store(:,j))
                nJ = temp_store_nI(:,j)

                ! If on the diagonal of the Hamiltonian.
                if (j == i + determ_displs(iProcIndex)) then
                    hamiltonian_row(j) = hphf_diag_helement(nI, SpawnedParts(:,i)) - Hii
                    core_ham_diag(i) = hamiltonian_row(j)
                    ! We calculate and store the diagonal matrix element at
                    ! this point for later access.
                    if (.not. tReadPops) &
                        call set_det_diagH(i, Real(hamiltonian_row(j), dp))
                    ! Always include the diagonal elements.
                    row_size = row_size + 1
                else
                    tmp = ieor(SpawnedParts(0:NIfD,i), temp_store(0:NIfD,j))
                    tmp = iand(SpawnedParts(0:NIfD,i), tmp)
                    IC = CountBits(tmp, NIfD)

                    if ( IC <= maxExcit .or. ((.not. CS_I) .and. (.not. cs(j))) ) then
                        hamiltonian_row(j) = hphf_off_diag_helement_opt(nI, SpawnedParts(:,i), temp_store(:,j), IC, CS_I, cs(j))
                        if (abs(hamiltonian_row(j)) > 0.0_dp) row_size = row_size + 1
                    end if
                end if

            end do

            ! Now we know the number of non-zero elements in this row of the Hamiltonian, so allocate it.
            call allocate_sparse_ham_row(sparse_core_ham, i, row_size, "sparse_core_ham", SparseCoreHamilTags(:,i))

            sparse_core_ham(i)%elements = 0.0_dp
            sparse_core_ham(i)%positions = 0
            sparse_core_ham(i)%num_elements = row_size

            counter = 1
            do j = 1, determ_space_size
                ! If non-zero or a diagonal element.
                if (abs(hamiltonian_row(j)) > 0.0_dp .or. (j == i + determ_displs(iProcIndex)) ) then
                    sparse_core_ham(i)%positions(counter) = j
                    sparse_core_ham(i)%elements(counter) = hamiltonian_row(j)
                    counter = counter + 1
                end if
                if (counter == row_size + 1) exit
            end do

        end do

        ! Don't time the mpi_barrier call, because if we did then we wouldn't
        ! be able separate out some of the core Hamiltonian creation time from
        ! the MPIBarrier calls in the main loop.
        call MPIBarrier(ierr, tTimeIn=.false.)

        deallocate(temp_store, stat=ierr)
        call LogMemDealloc(t_r, TempStoreTag, ierr)
        deallocate(hamiltonian_row, stat=ierr)
        call LogMemDealloc(t_r, HRTag, ierr)
        deallocate(temp_store_nI, stat=ierr)
        deallocate(cs, stat=ierr)

    end subroutine calc_determ_hamil_sparse_hphf

    subroutine calc_approx_hamil_sparse_hphf()

        use FciMCData, only: core_space

        integer :: i, j, row_size, counter, ierr
        integer :: nI(nel), nJ(nel)
        integer, allocatable :: temp_store_nI(:,:)
        HElement_t(dp), allocatable, dimension(:) :: hamiltonian_row
        character(len=*), parameter :: t_r = "calc_approx_hamil_sparse_hphf"

        integer(n_int) :: tmp(0:NIfD)
        integer :: IC
        logical :: CS_I, var_state_i, var_state_j
        logical, allocatable :: cs(:)

        allocate(approx_ham(determ_sizes(iProcIndex)), stat=ierr)
        allocate(hamiltonian_row(determ_space_size), stat=ierr)
        allocate(temp_store_nI(nel, determ_space_size), stat=ierr)
        allocate(cs(determ_space_size), stat=ierr)

        do i = 1, determ_space_size
            call decode_bit_det(temp_store_nI(:,i), core_space(:,i))
            cs(i) = TestClosedShellDet(core_space(:,i))
        end do

        ! Loop over all deterministic states on this processor.
        do i = 1, determ_sizes(iProcIndex)

            !call decode_bit_det(nI, SpawnedParts(:, i))
            nI = temp_store_nI(:, i + determ_displs(iProcIndex))

            var_state_i = is_var_state(core_space(:,i + determ_displs(iProcIndex)), nI)

            row_size = 0
            hamiltonian_row = 0.0_dp

            CS_I = cs(i + determ_displs(iProcIndex))

            ! Loop over all deterministic states.
            do j = 1, determ_space_size

                !call decode_bit_det(nJ, core_space(:,j))
                nJ = temp_store_nI(:,j)

                ! If on the diagonal of the Hamiltonian.
                if (j == i + determ_displs(iProcIndex)) then
                    hamiltonian_row(j) = hphf_diag_helement(nI, core_space(:, i + determ_displs(iProcIndex))) - Hii
                    !core_ham_diag(i) = hamiltonian_row(j)
                    ! We calculate and store the diagonal matrix element at
                    ! this point for later access.
                    if (.not. tReadPops) &
                        call set_det_diagH(i, Real(hamiltonian_row(j), dp))
                    ! Always include the diagonal elements.
                    row_size = row_size + 1
                else
                    var_state_j = is_var_state(core_space(:,j), nJ)

                    ! Only add a matrix element if both states are variational
                    ! states, or if one of them is (the rectangular portion of
                    ! H connected var_space to the space of connections to it).
                    if (var_state_i .or. var_state_j) then
                        tmp = ieor(core_space(0:NIfD,i+determ_displs(iProcIndex)), core_space(0:NIfD,j))
                        tmp = iand(core_space(0:NIfD,i+determ_displs(iProcIndex)), tmp)
                        IC = CountBits(tmp, NIfD)

                        if ( IC <= maxExcit .or. ((.not. CS_I) .and. (.not. cs(j))) ) then
                            hamiltonian_row(j) = hphf_off_diag_helement_opt(nI, core_space(:,i+determ_displs(iProcIndex)), &
                                                                             core_space(:,j), IC, CS_I, cs(j))

                            if (abs(hamiltonian_row(j)) > 0.0_dp) row_size = row_size + 1
                        end if
                    end if
                end if

            end do

            ! Now we know the number of non-zero elements in this row of the Hamiltonian, so allocate it.
            allocate(approx_ham(i)%elements(row_size), stat=ierr)
            allocate(approx_ham(i)%positions(row_size), stat=ierr)

            approx_ham(i)%elements = 0.0_dp
            approx_ham(i)%positions = 0
            approx_ham(i)%num_elements = row_size

            counter = 1
            do j = 1, determ_space_size
                ! If non-zero or a diagonal element.
                if (abs(hamiltonian_row(j)) > 0.0_dp .or. (j == i + determ_displs(iProcIndex)) ) then
                    approx_ham(i)%positions(counter) = j
                    approx_ham(i)%elements(counter) = hamiltonian_row(j)
                    counter = counter + 1
                end if
                if (counter == row_size + 1) exit
            end do

        end do

        ! Don't time the mpi_barrier call, because if we did then we wouldn't
        ! be able separate out some of the core Hamiltonian creation time from
        ! the MPIBarrier calls in the main loop.
        call MPIBarrier(ierr, tTimeIn=.false.)

        deallocate(hamiltonian_row, stat=ierr)
        deallocate(temp_store_nI, stat=ierr)
        deallocate(cs, stat=ierr)

    end subroutine calc_approx_hamil_sparse_hphf

    subroutine allocate_sparse_ham_row(sparse_matrix, row, sparse_row_size, sparse_matrix_name, sparse_tags)

        ! Allocate a single row and add it to the memory manager.

        type(sparse_matrix_real), intent(inout) :: sparse_matrix(:)
        integer, intent(in) :: row, sparse_row_size
        character(len=*), intent(in) :: sparse_matrix_name
        integer(TagIntType), intent(inout) :: sparse_tags(2)
        integer :: ierr
        character(len=1024) :: string_row
        character(len=1024) :: var_name
        character(len=*), parameter :: t_r = "allocate_sparse_ham_row"

        write (string_row, '(I10)') row

        var_name = trim(sparse_matrix_name)//"_"//trim(string_row)//"_elements"
        allocate(sparse_matrix(row)%elements(sparse_row_size), stat=ierr)
        call LogMemAlloc(var_name, sparse_row_size, 8, t_r, sparse_tags(1), ierr)

        var_name = trim(sparse_matrix_name)//"_"//trim(string_row)//"_positions"
        allocate(sparse_matrix(row)%positions(sparse_row_size), stat=ierr)
        call LogMemAlloc(var_name, sparse_row_size, bytes_int, t_r, sparse_tags(2), ierr)

    end subroutine allocate_sparse_ham_row

    subroutine deallocate_sparse_ham(sparse_matrix, sparse_tags)

        ! Deallocate the whole array, and remove all rows from the memory manager.

        type(sparse_matrix_real), intent(inout), allocatable :: sparse_matrix(:)
        integer(TagIntType), intent(inout), allocatable :: sparse_tags(:,:)
        integer :: sparse_matrix_size, i, ierr
        character(len=*), parameter :: t_r = "deallocate_sparse_ham"

        sparse_matrix_size = size(sparse_matrix)

        do i = sparse_matrix_size, 1, -1

            deallocate(sparse_matrix(i)%elements, stat=ierr)
            !call LogMemDealloc(t_r, sparse_tags(1,i), ierr)

            deallocate(sparse_matrix(i)%positions, stat=ierr)
            !call LogMemDealloc(t_r, sparse_tags(2,i), ierr)

        end do

        if (allocated(sparse_tags)) deallocate(sparse_tags)
        if (allocated(sparse_matrix)) deallocate(sparse_matrix)

    end subroutine deallocate_sparse_ham

    subroutine deallocate_sparse_matrix_int(sparse_mat)

        type(sparse_matrix_int), intent(inout), allocatable :: sparse_mat(:)

        integer :: i, ierr

        if (allocated(sparse_mat)) then
            do i = 1, size(sparse_mat)
                if (allocated(sparse_mat(i)%elements)) then
                    deallocate(sparse_mat(i)%elements, stat=ierr)
                    if (ierr /= 0) write(6,'("Error when deallocating sparse matrix elements array:",1X,i8)') ierr
                end if
                if (allocated(sparse_mat(i)%positions)) then
                    deallocate(sparse_mat(i)%positions, stat=ierr)
                    if (ierr /= 0) write(6,'("Error when deallocating sparse matrix positions array:",1X,i8)') ierr
                end if
            end do

            deallocate(sparse_mat, stat=ierr)
            if (ierr /= 0) write(6,'("Error when deallocating sparse matrix array:",1X,i8)') ierr
        end if

    end subroutine deallocate_sparse_matrix_int

    subroutine deallocate_core_hashtable(ht)

        type(core_hashtable), intent(inout), allocatable :: ht(:)

        integer :: i, ierr

        if (allocated(ht)) then
            do i = 1, size(ht)
                if (allocated(ht(i)%ind)) then
                    deallocate(ht(i)%ind, stat=ierr)
                    if (ierr /= 0) write(6,'("Error when deallocating core hashtable ind array:",1X,i8)') ierr
                end if
            end do

            deallocate(ht, stat=ierr)
            if (ierr /= 0) write(6,'("Error when deallocating core hashtable:",1X,i8)') ierr
        end if

    end subroutine deallocate_core_hashtable

    subroutine deallocate_trial_hashtable(ht)

        type(trial_hashtable), intent(inout), allocatable :: ht(:)

        integer :: i, ierr

        if (allocated(ht)) then
            do i = 1, size(ht)
                if (allocated(ht(i)%states)) then
                    deallocate(ht(i)%states, stat=ierr)
                    if (ierr /= 0) write(6,'("Error when deallocating trial hashtable states array:",1X,i8)') ierr
                end if
            end do

            deallocate(ht, stat=ierr)
            if (ierr /= 0) write(6,'("Error when deallocating core hashtable:",1X,i8)') ierr
        end if

    end subroutine deallocate_trial_hashtable

    function is_var_state(ilut, nI) result (var_state)

        use FciMCData, only: var_space, var_space_size_int
        use hash, only: FindWalkerHash

        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer, intent(in) :: nI(:)
        integer(int64) :: hash_val
        integer(int64) :: pos
        logical :: var_state

        hash_val = FindWalkerHash(nI, var_space_size_int)

        call var_ht%callback_lookup(hash_val, pos, var_state, loc_verify)
    contains

        function loc_verify(ind) result(match)
            integer(int64), intent(in) :: ind
            logical :: match

            match = all(ilut(0:NIfDBO) == var_space(0:NIfDBO, ind) )

        end function loc_verify
    end function is_var_state

end module sparse_arrays
