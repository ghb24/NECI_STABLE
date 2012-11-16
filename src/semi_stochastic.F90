#include "macros.h"

module semi_stochastic

    use bit_rep_data, only: flag_deterministic, nIfDBO, nOffY, nIfY, nOffFlag, &
                            deterministic_mask
    use bit_reps, only: NIfD, NIfTot, decode_bit_det, encode_bit_rep
    use CalcData, only: tRegenDiagHEls
    use csf, only: csf_get_yamas, get_num_csfs, get_csf_bit_yama, csf_apply_yama
    use csf_data, only: iscsf, csf_orbital_mask
    use constants
    use DetBitOps, only: ilut_lt, ilut_gt, count_open_orbs
    use Determinants, only: get_helement
    use enumerate_excitations
    use FciMCData, only: HFDet, ilutHF, iHFProc, Hii, CurrentDets, CurrentH, &
                         deterministic_proc_sizes, deterministic_proc_indices, &
                         det_vector, result_det_vector, core_hamiltonian
    use hash , only : DetermineDetNode
    use hphf_integrals, only: hphf_diag_helement
    use Parallel_neci, only: iProcIndex, nProcessors, MPIBCast, MPIBarrier
    use sort_mod, only: sort
    use SystemData, only: nel, tHPHF, tCSFCore, tDeterminantCore, Stot, lms

    implicit none

contains

    subroutine init_semi_stochastic()

        ! Initialise the semi-stochastic information. This includes enumerating a list
        ! of all determinants or CSFs in the deterministic space and calculating and
        ! storing the resulting Hamiltonian matrix elements. The lists which will store
        ! the psip vectors in the deterministic space are also allocated.

        integer :: i, j, IC
        integer :: det_space_size

        ! Initialise the deterministic mask.
        deterministic_mask = 0
        deterministic_mask = ibset(deterministic_mask, flag_deterministic)

        allocate(deterministic_proc_sizes(0:nProcessors-1))
        allocate(deterministic_proc_indices(0:nProcessors-1))
        deterministic_proc_sizes = 0
        deterministic_proc_indices = 0

        ! Count the Hartree-Fock determinant (The Hartree-Fock determinant has already been
        ! added to the main list as usual in the subroutine InitFCIMCCalcPar).
        deterministic_proc_sizes(iHFProc) = 1

        ! The following subroutines call the enumerating subroutines to create all excitations
        ! and add these states to the main list, CurrentDets, on the correct processor. As
        ! they do this, they count the size of the deterministic space on each processor.
        if (tDeterminantCore) then
            call generate_determinants()
        else if (tCSFCore) then
            call generate_csfs()
        else
            call stop_all("init_semi_stochastic", "The nature of the core basis functions &
                             &has not been specified.")
        end if

        ! We now know the size of the deterministic space on each processor, so now find
        ! the total size of the space and also allocate vectors to store psip amplitudes
        ! and the deterministic Hamiltonian.
        det_space_size = sum(deterministic_proc_sizes)
        allocate(det_vector(deterministic_proc_sizes(iProcIndex)))
        allocate(result_det_vector(det_space_size))
        allocate(core_hamiltonian(deterministic_proc_sizes(iProcIndex), det_space_size))
        det_vector = 0.0_dp
        result_det_vector = 0.0_dp

        ! Calculate the indices in the full vector at which the various processors take
        ! over, relative to the first index position in the vector (ie, the first value
        ! in this vector will be 0, *not* 1 - this is required for mpiallgatherv later).
        do i = 0, nProcessors-1
            if (i == 0) then
                deterministic_proc_indices(i) = 0
            else
                deterministic_proc_indices(i) = deterministic_proc_indices(i-1) + &
                                                  deterministic_proc_sizes(i-1)
            end if
        end do

        ! Sort the list of basis states so that it is correctly ordered in the predefined
        ! order which is always kept throughout the simulation.
        call sort(CurrentDets(:,1:deterministic_proc_sizes(iProcIndex)), ilut_lt, ilut_gt)

        ! Calculate and store the deterministic Hamiltonian.
        call calculate_det_hamiltonian_normal()

        call stop_all("here","here")

    end subroutine init_semi_stochastic

    subroutine calculate_det_hamiltonian_normal()

        integer :: i, j, iproc, col_index, ierr
        integer :: nI(nel), nJ(nel)
        integer(n_int), allocatable, dimension(:,:) :: temp_store

        ! temp_store is storage space for bitstrings so that the Hamiltonian matrix
        ! elements can be calculated.
        allocate(temp_store(0:NIfTot, maxval(deterministic_proc_sizes)))

        ! Calcuation of the Hamiltonian matrix elements to be stored. Loop over each
        ! processor in order and broadcast the sorted vector of bitstrings from the
        ! processor to all other processors so that all matrix elements can be found.
        do iproc = 0, nProcessors-1

            ! If we are broadcasting from this processor next, transfer the bitstrings
            ! to the array temp_store.
            if (iproc == iProcIndex) temp_store(:,1:deterministic_proc_sizes(iproc)) = &
                                      CurrentDets(:,1:deterministic_proc_sizes(iproc))

            ! Perform the broadcasting to other all other processors.
            call MPIBCast(temp_store, size(temp_store), iproc)

            ! The index of the column before the first column of the block of the
            ! Hamiltonian currently being calculated.
            col_index = deterministic_proc_indices(iproc)

            ! Loop over all the elements in the block of the Hamiltonian corresponding
            ! to these two prcoessors.
            do i = 1, deterministic_proc_sizes(iProcIndex)
                do j = 1, deterministic_proc_sizes(iproc)

                    ! Get nI form of the basis functions.
                    call decode_bit_det(nI, CurrentDets(:, i))
                    call decode_bit_det(nJ, temp_store(:, j))

                    ! If on the diagonal of the Hamiltonian.
                    if ((iprocIndex == iproc) .and. (i == j)) then
                        core_hamiltonian(i, col_index + j) = &
                                             get_helement(nI, nJ, 0)
                        ! We calculate and store CurrentH at this point for ease.
                        if (.not.tRegenDiagHEls) CurrentH(1,i) = &
                                             core_hamiltonian(i, col_index + j) - Hii
                    else
                        core_hamiltonian(i, col_index + j) = &
                                         get_helement(nI, nJ, CurrentDets(:, i), &
                                                            temp_store(:, j))
                    end if

                end do
            end do

        end do

        call MPIBarrier(ierr)

        deallocate(temp_store)

    end subroutine calculate_det_hamiltonian_normal

    subroutine generate_determinants()

        type(excit_store), target :: gen_store
        integer(n_int) :: ilut(0:NIfTot)
        integer :: nI(nel)

        ! This condition tells the enumerating subroutines to initialise the loop.
        ilut(0) = -1

        ! Find the first single excitation.
        call enumerate_all_single_excitations (ilutHF, HFDet, ilut, gen_store)
        ! Subroutine to add this state to the spawned list if on this processor.
        call add_basis_state_to_list(ilut)

        ! When no more basis functions are found, this value is returned and the loop
        ! is exited.
        do while(ilut(0) /= -1)
            call enumerate_all_single_excitations (ilutHF, HFDet, ilut, gen_store)

            ! If a determinant is returned (if we did not find the final one last time.)
            if (ilut(0) /= -1) call add_basis_state_to_list(ilut)
        end do

        ! Now generate the double excitations...

        ilut(0) = -1
        ! The first double excitation.
        call enumerate_all_double_excitations (ilutHF, HFDet, ilut, gen_store)
        call add_basis_state_to_list(ilut)

        do while(ilut(0) /= -1)
            call enumerate_all_double_excitations (ilutHF, HFDet, ilut, gen_store)

            if (ilut(0) /= -1) call add_basis_state_to_list(ilut)
        end do

    end subroutine generate_determinants

    subroutine generate_csfs()

        type(excit_store), target :: gen_store
        integer(n_int) :: ilutHF_loc(0:NIfTot), ilut(0:NIfTot)
        integer :: HFDet_loc(nel), nI(nel)
        integer :: exflag, n_open, i
        integer :: num_csfs, max_num_csfs
        integer, allocatable, dimension(:,:) :: yama_symbols
        logical :: first_loop

        ! Setting the first two bits of this flag tells the generating subroutine to
        ! generate both single and double spatial excitations.
        exflag = 0
        exflag = ibset(exflag, 0)
        exflag = ibset(exflag, 1)

        ! For Stot /= 0, the HF state will be a CSF. For the purpose of
        ! generating all spatial orbitals, we just want a determinant, so use a
        ! state without the CSF information.
        HFdet_loc = iand(HFDet, csf_orbital_mask)
        ilutHF_loc = ilutHF
        ilutHF_loc(NOffY) = 0

        ! Can't possibly have more open orbitals than nel.
        ! TODO: Check that num_csfs *always* increases with n_open.
        max_num_csfs = get_num_csfs(nel, Stot)

        ! Use max_num_csfs to allocate the array of Yamanouchi symbols to be large enough.
        allocate(yama_symbols(max_num_csfs, nel))
        yama_symbols = 0

        ! Ensure ilut(0) \= -1 so that the loop can be entered.
        ilut(0) = 0
        first_loop = .true.
        do while (ilut(0) /= -1)

            if (first_loop) then
                ilut(0) = -1
                first_loop = .false.
            end if

            ! Generate the next spatial excitation (orbital configuration).
            call enumerate_spatial_excitations(ilutHF_loc, HFDet_loc, ilut, exflag, gen_store)

            if (ilut(0) == -1) exit

            ! Find the nI representation.
            call decode_bit_det(nI, ilut)

            ! Find the number of open orbitals.
            n_open = count_open_orbs(ilut)
            ! Find the number of CSFs (and hence Yamanouchi symbols) with this value of Stot for
            ! this orbital configuration.
            num_csfs = get_num_csfs(n_open, Stot)

            ! Enumerate the list of all possible Yamanouchi symbols for this orbital
            ! configuration.
            call csf_get_yamas(n_open, Stot, yama_symbols(1:num_csfs,1:n_open), num_csfs)

            ! Loop over all Yamanouchi symbols.
            do i = 1, num_csfs
                ! If n_open = 0 then we just have a determinant.
                if (n_open > 0) then
                    ! Encode the Yamanouchi symbol in nI representation.
                    call csf_apply_yama(nI, yama_symbols(i,1:n_open))
                    ! Encode the Yamanouchi symbol in the ilut representation.
                    call get_csf_bit_yama(nI, ilut(nOffY:nOffY+nIfY-1))
                end if
                ! Finally add the CSF to the CurrentDets list.
                call add_basis_state_to_list(ilut, nI)
            end do

        end do

    end subroutine generate_csfs

    subroutine add_basis_state_to_list(ilut, nI_in)

        ! This subroutine, called from init_semi_stochastic, takes a bitstring, finds
        ! the nI representation, decides if the basis state lives on this processor
        ! and, if so, adds it to the spawned list. It also increases the vector sizes
        ! being calculated.

        ! In: ilut - The determinant in a bitstring form.
        ! In (optional) : nI_in - A list of the occupied orbitals in the determinant.

        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer, optional :: nI_in(nel)
        integer :: nI(nel)
        integer(n_int) :: flags
        integer :: proc, sgn(lenof_sign)

        sgn = 0
        ! Flag to specify that these basis states are in the deterministic space.
        flags = 0
        flags = ibset(flags,flag_deterministic)

        ! Find the nI representation of determinant.
        if (present(nI_in)) then
            nI = nI_in
        else
            call decode_bit_det(nI, ilut)
        end if

        ! Find the processor which this state belongs to.
        proc = DetermineDetNode(nI,0)

        ! Increase the size of the deterministic space on the correct processor.
        deterministic_proc_sizes(proc) = deterministic_proc_sizes(proc) + 1

        ! If this determinant belongs to this processor, add it to the main list.
        if (proc == iProcIndex) call encode_bit_rep(CurrentDets(:, &
                       deterministic_proc_sizes(iProcIndex)), ilut(0:nIfDBO), sgn, flags)

    end subroutine add_basis_state_to_list

end module semi_stochastic
