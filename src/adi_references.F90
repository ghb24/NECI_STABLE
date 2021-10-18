#include "macros.h"
module adi_references
    use Parallel_neci
    use FciMCData, only: ilutRef, TotWalkers, CurrentDets, core_run
    use adi_data, only: ilutRefAdi, nRefs, nIRef, signsRef, &
                        tAdiActive, tSetupSIs, NoTypeN, tSetupSIs, &
                        tReferenceChanged, SIThreshold, tUseCaches, nIRef, signsRef, exLvlRef, tSuppressSIOutput, &
                        targetRefPop, lastAllNoatHF, lastNRefs, tVariableNRef, maxNRefs, minSIConnect, &
                        nIncoherentDets, nConnection, tWeightedConnections, tSignedRepAv
    use CalcData, only: InitiatorWalkNo
    use bit_rep_data, only: niftot, nifd, extract_sign
    use bit_reps, only: decode_bit_det
    use DetBitOps, only: FindBitExcitLevel, sign_gt, sign_lt
    use sort_mod, only: sort
    use constants
    use SystemData, only: nel, t_3_body_excits, tGUGA
    use util_mod, only: operator(.isclose.)

    implicit none
    private
    public :: initialize_c_caches, update_coherence_check, eval_coherence, &
        check_sign_coherence, update_first_reference, enable_adi, clean_adi, &
        adjust_nRefs, update_ref_signs, resize_ilutRefAdi, &
        print_reference_notification, nRefs, setup_reference_space, &
        reallocate_ilutRefAdi, reset_coherence_counter, update_reference_space


contains

    subroutine setup_reference_space(tPopPresent)
        implicit none
        logical, intent(in) :: tPopPresent
        if (tAdiActive) then
            call update_reference_space(tPopPresent)
        end if
    end subroutine setup_reference_space

!------------------------------------------------------------------------------------------!

    subroutine update_reference_space(tPopPresent)
        use adi_data, only: tReadRefs
        use LoggingData, only: ref_filename
        implicit none
        logical, intent(in) :: tPopPresent
        integer :: nRead
        logical :: tGen

        if (tAdiActive) then
            ! If we actually did something
            tGen = .false.
            ! First, generate the reference space of nRefs determinants from the population (if present)
            if (tPopPresent) then
                call generate_ref_space()
                tGen = .true.
                if (.not. tReadRefs) &
                    call print_reference_notification(1, nRefs, "Superinitiators from population")
            end if

            if (tReadRefs) then
                ! Then, add the references from the file
                call read_in_refs(ref_filename, nRead)
                ! If we added references, note this (this means that the nRefs was not
                ! specified and is taken from the file)
                if (nRead > nRefs) nRefs = nRead
                tGen = .true.

                ! Print the read-in references
                call print_reference_notification(1, nRead, "Read in superinitiators")
                call print_reference_notification(nRead + 1, nRefs, "Superinitiators from population")
            end if
            ! Then, add the product references
            if (.not. tGen) then
                ! If we did not do anything, only take one reference
                call reallocate_ilutRefAdi(1)
                ilutRefAdi(:, 1) = ilutRef(:, 1)
                nRefs = 1
                call print_reference_notification(1, 1, &
                                                  "Using only the reference determinant as superinitiator")
            end if

            call fill_adi_caches()

            if (nRefs > 0) tSetupSIs = .true.
            tReferenceChanged = .true.
            call reset_coherence_counter()
        end if

    end subroutine update_reference_space

!------------------------------------------------------------------------------------------!

    subroutine fixed_number_SI_generation()
        use semi_stoch_gen, only: generate_space_most_populated
        use semi_stoch_procs, only: GLOBAL_RUN
        use LoggingData, only: ref_filename, tWriteRefs
        implicit none
        integer(MPIArg) :: mpi_refs_found
        integer :: ierr, i, all_refs_found, refs_found
        integer(MPIArg) :: refs_found_per_proc(0:nProcessors - 1), refs_displs(0:nProcessors - 1)
        integer(n_int) :: ref_buf(0:NIfTot, maxNRefs)

        ! we need to be sure ilutRefAdi has the right size
        call reallocate_ilutRefAdi(maxNRefs)
        if (maxNRefs > 0) then
            ! Get the nRefs most populated determinants
            refs_found = 0
            call generate_space_most_populated(maxNRefs, .false., 1, ref_buf, refs_found, &
                                               GLOBAL_RUN, CurrentDets(0:NifTot, :), TotWalkers)
            ! Communicate the refs_found info
            mpi_refs_found = int(refs_found, MPIArg)
            call MPIAllGather(mpi_refs_found, refs_found_per_proc, ierr)
            all_refs_found = sum(refs_found_per_proc)
            if (all_refs_found /= maxNRefs) then
                write(stdout, *) "Warning: Less than ", maxNRefs, &
                    " superinitiators found, using only ", all_refs_found, " superinitiators"
            end if
            ! Set the number of SIs to the number actually found
            nRefs = all_refs_found
            ! Prepare communication of SIs
            refs_displs(0) = 0
            do i = 1, nProcessors - 1
                refs_displs(i) = refs_displs(i - 1) + refs_found_per_proc(i - 1)
            end do
            ! Store them on all processors
            call MPIAllGatherV(ref_buf(0:NIfTot, 1:refs_found), ilutRefAdi, &
                               refs_found_per_proc, refs_displs)

            if (tWriteRefs) call output_reference_space(ref_filename)
        else
            nRefs = 0
        end if
    end subroutine fixed_number_SI_generation

!------------------------------------------------------------------------------------------!

    subroutine read_in_refs(filename, nRead)
        use util_mod, only: get_free_unit
        implicit none
        integer, intent(out) :: nRead
        character(255), intent(in) :: filename
        integer :: iunit, i, stat
        integer(n_int) :: tmp(0:NIfTot)
        logical :: exists
        character(*), parameter :: this_routine = "read_in_refs"

        ! All procs read in the file, as all need all references
        inquire (file=trim(filename), exist=exists)
        if (.not. exists) call stop_all(this_routine, "No "//trim(filename)//" file detected.")
        iunit = get_free_unit()

        ! If nRefs == 1, we set it to the number of references in the file
        if (nRefs == 1) then
            ! Therefore, we scan the file once
            open(iunit, file=trim(filename), status='old')
            nRefs = 0
            do
                read(iunit, *, iostat=stat) tmp
                if (stat < 0) exit
                nRefs = nRefs + 1
            end do
            close(iunit)

            ! And allocate the ilutRefAdi accordingly
            call reallocate_ilutRefAdi(nRefs)
        end if

        open(iunit, file=trim(filename), status='old')

        ! Read in at most nRefs references
        nRead = 0
        do i = 1, nRefs
            ! Note that read-in references always have precedence over generated references
            read(iunit, *, iostat=stat) ilutRefAdi(:, i)
            ! If there are no more dets to be read, exit
            if (stat < 0) exit
            ! If we successfully read, log it
            nRead = nRead + 1
        end do

        close(iunit)
    end subroutine read_in_refs

!------------------------------------------------------------------------------------------!

    subroutine generate_ref_space()
        use LoggingData, only: ref_filename, tWriteRefs
        implicit none
        integer :: refs_found, all_refs_found
        integer(n_int) :: ref_buf(0:NIfTot, maxNRefs), si_buf(0:NIfTot, maxNRefs)

        if (NoTypeN > InitiatorWalkNo) then
            call get_threshold_based_SIs(ref_buf, refs_found)

            ! communicate the SIs
            call communicate_threshold_based_SIs(si_buf, ref_buf, refs_found, all_refs_found)

            ! And write the so-merged references to ilutRef
            ! we need to be sure ilutRefAdi has the right size
            nRefs = all_refs_found
            call reallocate_ilutRefAdi(nRefs)
            ilutRefAdi(0:NIfTot, 1:nRefs) = si_buf(0:NIfTot, 1:nRefs)
        else
            call fixed_number_SI_generation()
        end if

        if (iProcIndex == root) &
            write(stdout, *) "Getting superinitiators for all-doubs-initiators: ", nRefs, " SIs found"

        if (tWriteRefs) call output_reference_space(ref_filename)
    end subroutine generate_ref_space

!------------------------------------------------------------------------------------------!

    subroutine get_threshold_based_SIs(ref_buf, refs_found)
        implicit none
        integer(n_int), intent(out) :: ref_buf(0:NIfTot, maxNRefs)
        integer, intent(out) :: refs_found
        integer(n_int), allocatable :: tmp(:, :)

!       integer :: i, nBlocks
        integer(int64) :: i
        integer :: nBlocks

        integer, parameter :: blockSize = 5000
        real(dp) :: sgn(lenof_sign)
        real(dp) :: repAvSgn

        ref_buf = 0
        refs_found = 0
        nBlocks = 1
        repAvSgn = 0.0_dp

        allocate(tmp(0:NIfTot, blockSize))
        tmp = 0

        do i = 1, TotWalkers
            call extract_sign(CurrentDets(:, i), sgn)
            ! either compare the sum of the signed or unsigned walker numbers to the
            ! initiator threshold
            if (tSignedRepAv) then
                repAvSgn = abs(sum(sgn) / inum_runs)
            else
                repAvSgn = av_pop(sgn)
            end if
            if ((repAvSgn >= NoTypeN)) then
                refs_found = refs_found + 1

                ! If the temporary is full, resize it
                if (refs_found > nBlocks * blockSize) then
                    nBlocks = nBlocks + 1
                    call resize_ilut_list(tmp, (nBlocks - 1) * blockSize, nBlocks * blockSize)
                end if

                ! add the possible SI to the temporary
                tmp(:, refs_found) = CurrentDets(:, i)
            end if
        end do

        ! we only keep at most maxNRefs determinants
        if (refs_found > maxNRefs .and. NoTypeN > 1) then
            write(stdout, '(A,I5,A,I5,A,I5,A)') "On proc ", iProcIndex, " found ", refs_found, &
                " SIs, which is more than the maximum of ", maxNRefs, " - truncating"
            ! in case we found more, take the maxNRefs with the highest population
            call sort(tmp(0:NIfTot, 1:refs_found), sign_gt, sign_lt)
            ref_buf(:, 1:maxNRefs) = tmp(:, 1:maxNRefs)
            refs_found = maxNRefs
        else
            ! else, take all of them
            ref_buf(:, 1:refs_found) = tmp(:, 1:refs_found)
        end if

        deallocate(tmp)
    end subroutine get_threshold_based_SIs

!------------------------------------------------------------------------------------------!

    subroutine communicate_threshold_based_SIs(si_buf, ref_buf, refs_found, all_refs_found)
        use semi_stoch_procs, only: return_largest_indices
        implicit none
        integer(n_int), intent(out) :: si_buf(0:NIfTot, maxNRefs)
        integer(n_int), intent(in) :: ref_buf(0:NIfTot, maxNRefs)
        integer, intent(in) :: refs_found
        integer, intent(out) :: all_refs_found
        integer(n_int), allocatable :: mpi_buf(:, :)
        integer(MPIArg) :: refs_found_per_proc(0:nProcessors - 1), refs_displs(0:nProcessors - 1)
        integer(MPIArg) :: mpi_refs_found
        integer :: ierr, i
        integer :: largest_inds(maxNRefs)
        real(dp), allocatable :: buf_signs(:)
        real(dp) :: tmp_sgn(lenof_sign)

        si_buf = 0
        ! here, we gather the potential SIs found by the procs and gather them
        ! Communicate the refs_found info
        mpi_refs_found = int(refs_found, MPIArg)
        call MPIAllGather(mpi_refs_found, refs_found_per_proc, ierr)
        ! total number of SI candiates
        all_refs_found = sum(refs_found_per_proc)
        refs_displs(0) = 0
        do i = 1, nProcessors - 1
            refs_displs(i) = refs_displs(i - 1) + refs_found_per_proc(i - 1)
        end do
        ! Store them on all processors
        allocate(mpi_buf(0:NIfTot, all_refs_found), stat=ierr)
        call MPIAllGatherV(ref_buf(0:NIfTot, 1:refs_found), mpi_buf, refs_found_per_proc, refs_displs)

        ! now, if we have have more than the maximum number of potential SIs, take the
        ! most populated only
        if (all_refs_found > maxNRefs) then
            ! get the indices of the largest elements in the communicated buffer
            ! first, extract the signs
            allocate(buf_signs(all_refs_found), stat=ierr)
            do i = 1, all_refs_found
                call extract_sign(mpi_buf(:, i), tmp_sgn)
                buf_signs(i) = sum(abs(tmp_sgn))
            end do
            ! then get the indices
            call return_largest_indices(maxNRefs, all_refs_found, buf_signs, largest_inds)
            ! and then copy those elements to the output buffer
            do i = 1, maxNRefs
                si_buf(0:NIfTot, i) = mpi_buf(0:NIfTot, largest_inds(i))
            end do
            if (iProcIndex == root .and. NoTypeN > 1) &
                write(stdout, '(A,I5,A,I5,A)') "In total ", all_refs_found, &
                " SIs were found, which is more than the maximum number of ", &
                maxNRefs, " - truncating"
            ! make it look to the outside as though maxNRefs were found
            all_refs_found = maxNRefs
        else
            ! else, take all elements
            si_buf(0:NIfTot, 1:all_refs_found) = mpi_buf(0:NIfTot, 1:all_refs_found)
        end if

        deallocate(mpi_buf)
    end subroutine communicate_threshold_based_SIs

!------------------------------------------------------------------------------------------!

    subroutine output_reference_space(filename)
        use MPI_wrapper, only: root
        use util_mod, only: get_free_unit
        implicit none
        character(255), intent(in) :: filename
        integer :: iunit, i, ierr

        if (iProcIndex == root) then
            iunit = get_free_unit()
            open(iunit, file=filename, status='replace')
            ! write the references into the file
            do i = 1, nRefs
                write(iunit, *) ilutRefAdi(:, i)
            end do
            call neci_flush(iunit)
            close(iunit)
        end if
        call MPIBarrier(ierr)
    end subroutine output_reference_space

!------------------------------------------------------------------------------------------!

    subroutine print_reference_notification(iStart, iEnd, title, legend)
        implicit none
        integer, intent(in) :: iStart, iEnd
        character(*), intent(in) :: title
        logical, optional :: legend
        integer :: i

        if (iProcIndex == root .and. .not. tSuppressSIOutput) then
            ! print out the given SIs sorted according to population
            call sort(ilutRefAdi(0:NIfTot, iStart:iEnd), sign_gt, sign_lt)
            write(stdout, *) title
            if (present(legend)) write(stdout, "(4A25)") &
                ! TODO: Adapt legend for multiple runs
                "Determinant (bitwise)", "Excitation level", "Coherence parameter", "Number of walkers"
            do i = iStart, iEnd
                call print_bit_rep(ilutRefAdi(:, i))
            end do
        end if
    end subroutine print_reference_notification

!------------------------------------------------------------------------------------------!

    subroutine print_bit_rep(ilut)
        implicit none
        integer(n_int), intent(in) :: ilut(0:NIfTot)
        real(dp) :: temp_sgn(lenof_sign)
        integer :: j
        ! First, print the determinant (bitwise)
        call WriteDetBit(stdout, ilut, .false.)
        ! Then, the excitation level
        write(stdout, "(5X)", advance='no')
        write(stdout, "(G1.4)", advance='no') FindBitExcitLevel(ilut, ilutRef(:, 1))
        ! And the sign coherence parameter
        write(stdout, "(G16.7)", advance='no') get_sign_op(ilut)
        ! And then the sign
        call extract_sign(ilut, temp_sgn)
        do j = 1, lenof_sign
            write(stdout, "(G16.7)", advance='no') temp_sgn(j)
        end do

        write(stdout, '()')

    end subroutine print_bit_rep

!------------------------------------------------------------------------------------------!

    subroutine update_ref_signs()
        implicit none
        integer :: iRef

        do iRef = 1, nRefs
            call update_single_ref_sign(iRef)
        end do
    end subroutine update_ref_signs

!------------------------------------------------------------------------------------------!

    subroutine update_single_ref_sign(iRef)
        ! Get the signs of the references from the currentdets. Used both in the final
        ! output and in generating the product excitations
        use bit_reps, only: decode_bit_det, encode_sign
        use hash, only: hash_table_lookup
        use FciMCData, only: HashIndex
        use Parallel_neci, only: MPISumAll
        implicit none
        integer, intent(in) :: iRef
        integer:: nI(nel), hash_val, index
        real(dp) :: tmp_sgn(lenof_sign), mpi_sgn(lenof_sign)
        logical :: tSuccess

        call decode_bit_det(nI, ilutRefAdi(:, iRef))
        call hash_table_lookup(nI, ilutRefAdi(:, iRef), nifd, HashIndex, CurrentDets, &
                               index, hash_val, tSuccess)
        if (tSuccess) then
            call extract_sign(CurrentDets(:, index), tmp_sgn)
        else
            tmp_sgn = 0.0_dp
        end if
        ! This does the trick of communication: We need to get the sign to all
        ! processors, so just sum up the individual ones
        call MPISumAll(tmp_sgn, mpi_sgn)

        ! The sign contains information on all runs, so it is the same on all
        call encode_sign(ilutRefAdi(:, iRef), mpi_sgn)

    end subroutine update_single_ref_sign


    function check_sign_coherence(ilut, nI, ilut_sign, iRef, run) result(is_coherent)
        use Determinants, only: get_helement
        use SystemData, only: tHPHF
        use hphf_integrals, only: hphf_off_diag_helement
        use guga_excitations, only: calc_guga_matrix_element
        use guga_bitRepOps, only: CSF_Info_t
        use guga_data, only: ExcitationInformation_t
        implicit none
        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer, intent(in) :: iRef, run, nI(nel)
        real(dp), intent(in) :: ilut_sign(lenof_sign)
        logical :: is_coherent
        integer :: nJRef(nel)
        real(dp) :: signRef(lenof_sign)
#ifdef CMPLX_
        complex(dp) :: tmp
#endif
        type(ExcitationInformation_t) :: excitInfo
        HElement_t(dp) :: h_el

        is_coherent = .true.
        ! Obviously, we need the sign to compare coherence
        call extract_sign(ilutRefAdi(:, iRef), signRef)

        ! For getting the matrix element, we also need the determinants
        call decode_bit_det(nJRef, ilutRefAdi(:, iRef))

        ! Then, get the matrix element
        if (tHPHF) then
            h_el = hphf_off_diag_helement(nI, nJRef(:), ilut, ilutRefAdi(:, iRef))
        else if (tGUGA) then
            call calc_guga_matrix_element(ilut, CSF_Info_t(ilut), ilutRefAdi(:, iref), excitInfo, h_el, .true., 2)
        else
            h_el = get_helement(nI, nJRef(:), ilut, ilutRefAdi(:, iRef))
        end if
        ! if the determinants are not coupled, ignore them
        if (abs(h_el) < eps) return

        ! If the new ilut has sign 0, there is no need to do any check on this run
        if (mag_of_run(ilut_sign, run) > eps) then
#ifdef CMPLX_
            ! The complex coherence check is more effortive, so only do it in complex builds
            ! Get the relative phase of the signs
            tmp = cmplx(signRef(min_part_type(run)), signRef(max_part_type(run)), dp) / &
                  cmplx(ilut_sign(min_part_type(run)), ilut_sign(max_part_type(run)), dp)
            ! and compare it to that of H
            if (aimag(h_el * tmp) > eps .or. real(h_el * tmp) > 0.0_dp) then
                is_coherent = .false.
                return
            end if
#else
            ! For the real comparison, we just compare the signs
            if (signRef(run) * ilut_sign(run) * h_el > 0.0_dp) then
                is_coherent = .false.
                return
            end if
#endif
        end if
    end function check_sign_coherence

!------------------------------------------------------------------------------------------!

    subroutine initialize_c_caches(signedCache, unsignedCache, connections)
        implicit none
        HElement_t(dp), intent(out) :: signedCache
        real(dp), intent(out) :: unsignedCache
        integer, intent(out) :: connections

        ! We potentially want to add the diagonal terms here
        signedCache = 0.0_dp
        unsignedCache = 0.0_dp
        connections = 0
    end subroutine initialize_c_caches

!------------------------------------------------------------------------------------------!

    subroutine update_coherence_check(ilut, nI, i, signedCache, unsignedCache, connections)
        use SystemData, only: tHPHF
        use Determinants, only: get_helement
        use hphf_integrals, only: hphf_off_diag_helement
        use guga_excitations, only: calc_guga_matrix_element
        use guga_bitRepOps, only: CSF_Info_t
        use guga_data, only: ExcitationInformation_t
        implicit none
        integer, intent(in) :: nI(nel), i
        integer(n_int), intent(in) :: ilut(0:NIfTot)
        HElement_t(dp), intent(inout) :: signedCache
        real(dp), intent(inout) :: unsignedCache
        integer, intent(out) :: connections
        real(dp) :: i_sgn(lenof_sign)
        integer :: run
        HElement_t(dp) :: h_el, tmp
        type(ExcitationInformation_t) :: excitInfo

        ! TODO: Only if ilutRefAdi(:,i) is a SI on this run

        ! We require cached data here
        if (.not. tUseCaches) call fill_adi_caches
        ! First, get the matrix element
        if (tHPHF) then
            h_el = hphf_off_diag_helement(nI, nIRef(:, i), ilut, ilutRefAdi(:, i))
        else if (tGUGA) then
            call calc_guga_matrix_element(ilut, CSF_Info_t(ilut), ilutRefAdi(:, i), excitInfo, h_el, .true., 2)
        else
            h_el = get_helement(nI, nIRef(:, i), ilut, ilutRefAdi(:, i))
        end if
        ! Only proceed if the determinants are coupled
        if (abs(h_el) < eps) return

        ! Add tmp = Hij cj to the caches

        tmp = h_cast(0.0_dp)
        do run = 1, inum_runs
#ifdef CMPLX_
            tmp = tmp + h_el * cmplx(signsRef(min_part_type(run), i), &
                                     signsRef(max_part_type(run), i), dp)
#else
            tmp = h_el * signsRef(run, i)
#endif
        end do
        signedCache = signedCache + tmp
        unsignedCache = unsignedCache + abs(tmp)
        if (tWeightedConnections) then
            ! there is the option to have the connections weighted with
            ! the population
            i_sgn = signsRef(:, i)
            do run = 1, inum_runs
                connections = int(connections + mag_of_run(i_sgn, run) / inum_runs)
            end do
        else
            connections = connections + 1
        end if

    end subroutine update_coherence_check

!------------------------------------------------------------------------------------------!

    subroutine eval_coherence(signedCache, unsignedCache, sgn, connections, staticInit)
        ! Note that tweakcoherentdoubles and tavcoherentdoubles are mutually exclusive
        ! in the current implementation
        use adi_data, only: tWeakCoherentDoubles, tAvCoherentDoubles, coherenceThreshold, &
                            nConnection
        implicit none
        HElement_t(dp), intent(in) :: signedCache
        real(dp), intent(in) :: unsignedCache, sgn
        integer, intent(in) :: connections
        logical, intent(inout) :: staticInit

        ! Only need to check if we are looking at a double
        !if(unsignedCache > eps .and. connections>=minSIConnect) then
        if (connections < minSIConnect .or. &
            ! if the connections are weighted, we want to have at least
            ! minimum-number-of-connections*SI-threshold
            (tWeightedConnections .and. connections < minSIConnect * NoTypeN)) then
            staticInit = .false.
            if (unsignedCache > EPS) nConnection = nConnection + 1
        end if
        ! We disable superinitiator-related initiators if they fail the coherence check
        ! else, we leave it as it is
        if (tWeakCoherentDoubles) then
            if (abs(signedCache) < coherenceThreshold * unsignedCache) then
                staticInit = .false.
                nIncoherentDets = nIncoherentDets + 1
            end if
        end if

        ! If we do averaged coherence check, we check versus the sign of the ilut
        ! We recommend using both, av and weak
        if (tAvCoherentDoubles) then
            if (real(signedCache * sgn, dp) > 0.0_dp) then
                staticInit = .false.
                nIncoherentDets = nIncoherentDets + 1
            end if
        end if

    end subroutine eval_coherence

!------------------------------------------------------------------------------------------!

    function get_sign_op(ilut) result(xi)
        ! compute the sign problem indicator xi for a given determinant
        implicit none
        integer(n_int), intent(in) :: ilut(0:NIfTot)
        real(dp) :: xi
        integer :: iRef, nI(nel), exLevel
        real(dp) :: unsignedCache
        HElement_t(dp) :: signedCache
#ifdef DEBUG_
        character(*), parameter :: this_routine = "get_sign_op_run"
#endif
        integer :: connections
        ASSERT(.not. t_3_body_excits)

        call initialize_c_caches(signedCache, unsignedCache, connections)
        ! Sum up all Hij cj for all superinitiators j
        do iRef = 1, nRefs
            exLevel = FindBitExcitLevel(ilutRefAdi(:, iRef), ilut)
            ! Of course, only singles/doubles of ilut can contribute
            if (exLevel < 3) then
                call decode_bit_det(nI, ilut)
                call update_coherence_check(ilut, nI, iRef, &
                                            signedCache, unsignedCache, connections)
            end if
        end do
        ! Get the ratio of the caches, which is xi
        if (abs(unsignedCache) > eps .and. connections >= minSIConnect) then
            xi = abs(signedCache) / unsignedCache
        else
            xi = 0.0_dp
        end if
    end function get_sign_op

!------------------------------------------------------------------------------------------!

    subroutine adjust_nRefs()
        ! update the number of SIs used
        use adi_data, only: targetRefPopTol, tSingleSteps
        use FciMCData, only: AllNoatHF
        implicit none
        real(dp) :: cAllNoatHF
        integer :: nRefsOld

        if (tVariableNRef) then

            nRefsOld = nRefs
            ! average reference population
            cAllNoatHF = sum(AllNoatHF) / inum_runs
            ! check if we reached the destination
            if (abs(targetRefPop - cAllNoatHF) < targetRefPopTol) return

            if (tSingleSteps) then
                if (cAllNoatHF < targetRefPop) then
                    ! we cannot go below 0 SIs, if we hit 0, disable the adi option
                    if (nRefs > 1) then
                        nRefs = nRefs - 1
                    else
                        ! might change the behaviour here, to whatever makes sense
                        tAdiActive = .false.
                    end if
                end if
                if (cAllNoatHF > targetRefPop) nRefs = nRefs + 1
                tSingleSteps = .false.
            else
                ! do the extrapolationx
                nRefs = guess_target_nref()
            end if

            ! store the current values for the next step
            lastAllNoatHF = cAllNoatHF
            lastNRefs = nRefsOld

            write(stdout, *) "Now at ", nRefs, " superinitiators"

        end if
    end subroutine adjust_nRefs

!------------------------------------------------------------------------------------------!

    function guess_target_nref() result(target_nref)
        ! try to extrapolate how many SIs are required for a given target population
        use FciMCData, only: AllNoatHF
        implicit none
        integer :: target_nref
        real(dp) :: alpha, beta, cAllNoatHF

        target_nref = nRefs
        cAllNoatHF = sum(AllNoatHF) / inum_runs
        ! if the reference population did not change, nothing to do
        if (lastAllNoatHF.isclose.cAllNoatHF) return

        ! if we did not change the SI threshold, nothing to do
        if (lastNRefs == nRefs) return
        beta = log(lastAllNoatHF / cAllNoatHF) / (lastNRefs - nRefs)
        alpha = lastAllNoatHF * exp(log(lastAllNoatHF / cAllNoatHF) * (lastNRefs / (lastNRefs - nRefs)))

        ! here, we really assume exponential behaviour and extrapolate the nRefs that would then
        ! give the reference population closest to the target
        target_nref = int(log(targetRefPop / alpha) / beta)

    end function guess_target_nref

!------------------------------------------------------------------------------------------!

    subroutine fill_adi_caches()
        implicit none
        integer :: iRef

        call allocate_adi_caches()

        do iRef = 1, nRefs
            call decode_bit_det(nIRef(:, iRef), ilutRefAdi(:, iRef))
            call extract_sign(ilutRefAdi(:, iRef), signsRef(:, iRef))
            exLvlRef(iRef) = FindBitExcitLevel(ilutRef(:, 1), ilutRefAdi(:, iRef))
        end do

        tUseCaches = .true.
    end subroutine fill_adi_caches
!------------------------------------------------------------------------------------------!

    subroutine allocate_adi_caches()
        implicit none

        call deallocate_adi_caches()

        allocate(nIRef(nel, nRefs))
        allocate(signsRef(lenof_sign, nRefs))
        allocate(exLvlRef(nRefs))
    end subroutine allocate_adi_caches

!------------------------------------------------------------------------------------------!

    subroutine deallocate_adi_caches()
        implicit none

        if (allocated(nIRef)) deallocate(nIRef)
        if (allocated(signsRef)) deallocate(signsRef)
        if (allocated(exLvlRef)) deallocate(exLvlRef)
    end subroutine deallocate_adi_caches

!------------------------------------------------------------------------------------------!

    subroutine resize_ilutRefAdi(newSize)
        implicit none
        integer, intent(in) :: newSize

        ! see resize_ilut_list for the implementation
        call resize_ilut_list(ilutRefAdi, nRefs, newSize)
        ! The new number of references
        nRefs = newSize

        tUseCaches = .false.

    end subroutine resize_ilutRefAdi

!------------------------------------------------------------------------------------------!

    subroutine resize_ilut_list(list, oldSize, newSize)
        implicit none
        integer, intent(in) :: newSize, oldSize
        integer(n_int), allocatable, intent(inout) :: list(:, :)
        integer(n_int) :: tmp(0:NIfTot, oldSize)
        integer :: ierr
        logical :: tUseTmp
        character(*), parameter :: this_routine = "resize_ilut_list"

        ! We store the current list in a temporary, if existent
        tUseTmp = .false.
        ! If the list is already allocated, clear it
        if (allocated(list)) then
            tmp(:, :) = list(:, :)
            tUseTmp = .true.
            deallocate(list)
        end if

        allocate(list(0:NIfTot, newSize), stat=ierr)
        if (ierr /= 0) call stop_all(this_routine, "Not enough memory to resize ilut list")
        list = 0

        ! Write the old entries of list into the new memory
        if (oldSize < newSize) then
            list(:, 1:oldSize) = tmp(:, 1:oldSize)
        else
            list(:, 1:newSize) = tmp(:, 1:newSize)
        end if

    end subroutine resize_ilut_list

!------------------------------------------------------------------------------------------!

    subroutine reallocate_ilutRefAdi(size)
        implicit none
        integer, intent(in) :: size

        ! reallocate first the ilutref itself
        if (allocated(ilutRefAdi)) deallocate(ilutRefAdi)
        allocate(ilutRefAdi(0:NIfTot, size))
        ilutRefAdi = 0_n_int

        tSetupSIs = .false.
        tUseCaches = .false.

    end subroutine reallocate_ilutRefAdi

!------------------------------------------------------------------------------------------!


!------------------------------------------------------------------------------------------!

    subroutine enable_adi
        use adi_data, only: tAllDoubsInitiators, tAdiActive
        implicit none

        write(stdout, '()')
        write(stdout, *) "Setting all double excitations to initiators"
        write(stdout, *) "Using static references"
        write(stdout, *) "Further notification on additional references will be given"
        write(stdout, '()')
        tAllDoubsInitiators = .true.
        tAdiActive = .true.

        ! tSinglePartPhase = .false.

    end subroutine enable_adi

!------------------------------------------------------------------------------------------!

    subroutine reset_coherence_counter()
        use adi_data, only: nCoherentDoubles
        implicit none

        nCoherentDoubles = 0
        nIncoherentDets = 0
        nConnection = 0
    end subroutine reset_coherence_counter

!------------------------------------------------------------------------------------------!

    subroutine update_first_reference()
        use FciMCData, only: tSinglePartPhase
        implicit none

        call setup_reference_space(all(.not. tSinglePartPhase))
    end subroutine update_first_reference

!------------------------------------------------------------------------------------------!

    subroutine clean_adi()
        implicit none

        call deallocate_adi_caches()
        if (allocated(ilutRefAdi)) deallocate(ilutRefAdi)
    end subroutine clean_adi

end module adi_references
