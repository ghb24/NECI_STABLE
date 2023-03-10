#include "macros.h"

module enumerate_excitations

    use SystemData, only: tReltvy, t_k_space_hubbard, t_new_real_space_hubbard

    use util_mod, only: operator(.div.)

    use bit_rep_data, only: NIfD, NIfTot, nifguga

    use bit_reps, only: decode_bit_det

    use constants

    use util_mod, only: stop_all

    use DetBitOps, only: IsAllowedHPHF, TestClosedShellDet, EncodeBitDet

    use FciMCData, only: SpinInvBRR

    use HPHFRandExcitMod, only: FindExcitBitDetSym

    use Parallel_neci, only: MPISumAll

    use sort_mod, only: sort

    use SymData, only: nSymLabels, SymTable, SymLabels, SymClasses

    use SymExcit3, only: GenExcitations3

    use SymExcit4, only: GenExcitations4, ExcitGenSessionType

    use SymExcitDataMod

    use sym_general_mod

    use SystemData, only: nel, nBasis, G1, tFixLz, Arr, Brr, tHPHF, tHub, &
                          tUEG, tKPntSym, tReal, tUseBrillouin, tGUGA, tReltvy

    use guga_data, only: tag_excitations
    use guga_bitRepOps, only: CSF_Info_t

    use MemoryManager, only: LogMemDealloc

#ifndef CMPLX_
    use lattice_models_utils, only: gen_all_excits_k_space_hubbard, &
                                    gen_all_excits_r_space_hubbard
#endif

    better_implicit_none

    ! This type allows data in the enumerating subroutines to be stored
    ! such that when the subroutine is next called it can begin where
    ! it left of and go to the next excitation.
    type excit_store
        integer :: i, j
        integer :: orb1, orb2
        integer :: sym_ind1, sym_ind2
        integer :: ind1, ind2
        integer :: count1, count2
        integer :: sym1, sym2
        integer :: norb_sym
        integer :: sym_prod
        integer :: npairs
        logical :: gen_singles
    end type

    type simple_excit_store
        integer :: i, j
    end type

    type(ExcitGenSessionType) :: session

contains


    subroutine init_generate_connected_space(nI, ex_flag, tAllExcitFound, excit, excit_gen, nstore, tTempUseBrill)

        use SymExcit2, only: gensymexcitit2par_worker

        integer, intent(in) :: nI(nel)
        integer, intent(out) :: ex_flag
        logical, intent(out) :: tAllExcitFound
        integer, intent(out) :: excit(2, 2)
        integer, allocatable, intent(out) :: excit_gen(:)
        integer, intent(out) :: nStore(6)
        logical, intent(out) :: tTempUseBrill

        integer :: iMaxExcit, nExcitMemLen(1), nJ(nel), ierr
        character(*), parameter :: t_r = 'init_generate_connected_space'

        ! Which excitations levels should be considered?
        ! Singles only (1), doubles only (2) or singles and doubles (3).
        if (tHub) then
            if (tReal) then
                ex_flag = 1
            else
                ex_flag = 2
            end if
        else if (tUEG) then
            ex_flag = 2
        else
            ex_flag = 3
        end if

        tAllExcitFound = .false.

        if (tKPntSym) then
            ! We have to ensure that brillouins theorem isn't on for the
            ! excitation generator.
            if (tUseBrillouin) then
                tTempUseBrill = .true.
                tUseBrillouin = .false.
            else
                tTempUseBrill = .false.
            end if

            iMaxExcit = 0
            nStore = 0

            call GenSymExcitIt2Par_worker(nI, nel, G1, nBasis, .true., nExcitMemLen, nJ, &
                                          iMaxExcit, nStore, ex_flag, 1, nel)

            allocate(excit_gen(nExcitMemLen(1)), stat=ierr)
            if (ierr /= 0) call stop_all(t_r, "Problem allocating excitation generator.")
            excit_gen = 0

            call GenSymExcitIt2Par_worker(nI, nel, G1, nBasis, .true., excit_gen, nJ, &
                                          iMaxExcit, nStore, ex_flag, 1, nel)
        else
            excit = 0
        end if

    end subroutine init_generate_connected_space

    subroutine generate_connected_space(original_space_size, original_space, &
                                        connected_space_size, connected_space, tSinglesOnlyOpt)

        ! A wrapper function to call the correct generation routine.

        integer, intent(in) :: original_space_size
        integer(n_int), intent(in) :: original_space(0:, :)
        integer, intent(inout) :: connected_space_size
        integer(n_int), optional, intent(out) :: connected_space(0:, :)
        logical, intent(in), optional :: tSinglesOnlyOpt
        character(*), parameter :: this_routine = "generate_connected_space"

        if (tKPntSym) then
            call generate_connected_space_kpnt(original_space_size, original_space, &
                                               connected_space_size, connected_space, tSinglesOnlyOpt)
        else
            call generate_connected_space_normal(original_space_size, original_space, &
                                                 connected_space_size, connected_space, tSinglesOnlyOpt)
        end if

    end subroutine generate_connected_space

    subroutine generate_connected_space_normal(original_space_size, original_space, &
                                               connected_space_size, connected_space, tSinglesOnlyOpt)

        use Symexcit4, only: NewParentDet
        ! This routine either counts or generates all the determinants connected to those in
        ! original_space. If connected_space is not present then they will only be counted,
        ! else they will be stored in connected_space. If tSinglesOnlyOpt is present and
        ! also .true. then actually this routine will only generate all single excitations
        ! (regardless of the system being studied), otherwise it will generate all connected
        ! determinants.

        use guga_bitRepOps, only: convert_ilut_toGUGA, convert_ilut_toNECI
        use guga_excitations, only: actHamiltonian
        use SystemData, only: tGUGA
        integer :: nexcit, j
        integer(n_int), allocatable :: excitations(:, :)
        integer(n_int) :: ilutG(0:nifguga)

        integer, intent(in) :: original_space_size
        integer(n_int), intent(in) :: original_space(0:, :)
        integer, intent(inout) :: connected_space_size
        integer(n_int), optional, intent(out) :: connected_space(0:, :)
        logical, intent(in), optional :: tSinglesOnlyOpt
        character(*), parameter :: this_routine = "generate_connection_normal"

        integer(n_int) :: ilutJ(0:NIfTot)
        integer :: nI(nel), nJ(nel)
        integer :: i, excit(2, 2), ex_flag
        integer, allocatable :: excit_gen(:)
        integer :: nStore(6)
        logical :: tAllExcitFound, tStoreConnSpace, tSinglesOnly, tTempUseBrill
        integer :: n_excits
        integer(n_int), allocatable :: temp_dets(:, :)

        if (present(connected_space)) then
            tStoreConnSpace = .true.
        else
            tStoreConnSpace = .false.
        end if

        tSinglesOnly = .false.
        if (present(tSinglesOnlyOpt)) then
            if (tSinglesOnlyOpt) tSinglesOnly = .true.
        end if

        connected_space_size = 0

        ! Over all the states in the original list:
        do i = 1, original_space_size

            call decode_bit_det(nI, original_space(0:NIfTot, i))

            ! do the GUGA changes here, I want to do all the excitations from
            ! the currently looped over original_space(:,i)
            ! i think i still want to do this this way, since the dets
            ! implementation is really akward..
            if (tGUGA) then
                ! in GUGA don't do the tSinglesOnly option
                ASSERT(.not. tSinglesOnly)

                ! only STORE the excitations if the proper flag is set,
                ! otherwise only, increase the counter for the connected space
                ! why is this done??
                call convert_ilut_toGUGA(original_space(:, i), ilutG)

                call actHamiltonian(ilutG, CSF_Info_t(ilutG), excitations, nexcit)

                ! and if store flag is present:
                if (tStoreConnSpace) then
                    do j = 1, nexcit
                        call convert_ilut_toNECI(excitations(:, j), &
                                                 connected_space(:, connected_space_size + j))
                    end do
                end if

                ! update connected_space_size afterwards
                connected_space_size = connected_space_size + nexcit

                deallocate(excitations)
                call LogMemDealloc(this_routine, tag_excitations)

            else
                if (t_new_real_space_hubbard) then

#ifdef CMPLX_
                    call stop_all(this_routine, "not implemented for complex")
#else
                    call gen_all_excits_r_space_hubbard(nI, n_excits, temp_dets)
#endif

                    if (tStoreConnSpace) then
                        connected_space(0:nifd, connected_space_size + 1:connected_space_size + n_excits) &
                            = temp_dets(0:nifd, :)
                    end if

                    connected_space_size = connected_space_size + n_excits

                else

                    call NewParentDet(session)

                    call init_generate_connected_space(nI, ex_flag, tAllExcitFound, excit, excit_gen, nstore, tTempUseBrill)

                    if (tSinglesOnly) ex_flag = 1

                    do while (.true.)

                        call generate_connection_normal(nI, original_space(:, i), nJ, ilutJ, ex_flag, excit, &
                                                        tAllExcitFound, ncon=connected_space_size)
                        if (tAllExcitFound) exit

                        if (tStoreConnSpace) connected_space(0:NIfD, connected_space_size) = ilutJ(0:NIfD)

                    end do
                end if

            end if ! tGUGA
        end do

    end subroutine generate_connected_space_normal

    subroutine generate_connection_normal(nI, ilutI, nJ, ilutJ, ex_flag, excit, &
                                          tAllExcitFound, hel, ncon)

        use procedure_pointers, only: get_conn_helement
        use SymExcit4, only: GenExcitations4

        integer :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel)
        integer(n_int), intent(out) :: ilutJ(0:NIfTot)
        integer, intent(inout) :: ex_flag
        integer, intent(inout) :: excit(2, 2)
        logical, intent(inout) :: tAllExcitFound
        HElement_t(dp), optional, intent(out) :: hel
        integer, optional, intent(inout) :: ncon

        integer(n_int) :: ilut_tmp(0:NIfTot)
        integer :: ic
        logical :: tParity
        HElement_t(dp) :: hel_unused

        ! Generate the next determinant.
        if (tReltvy) then
            call GenExcitations4(session, nI, nJ, ex_flag, excit, tParity, tAllExcitFound, .false.)
        else
            call GenExcitations3(nI, ilutI, nJ, ex_flag, excit, tParity, &
                                 tAllExcitFound, .false.)
        end if
        if (tAllExcitFound) return

        ! Encode nJ in ilutJ.
        call EncodeBitDet(nJ, ilutJ)

        ! If using HPHFs, and if this isn't the allowed state for this HPHF, find
        ! the the allowed state and continue with this as ilut.
        if (tHPHF .and. (.not. TestClosedShellDet(ilutJ))) then
            if (.not. IsAllowedHPHF(ilutJ(0:NIfD))) then
                ilut_tmp = ilutJ
                call FindExcitBitDetSym(ilut_tmp, ilutJ)
            end if
        end if

        if (present(HEl)) then
            ! Map ex_flag values (1, 2 or 3) to 1, 2 and 1, respectively.
            ic = 2 - mod(ex_flag, 2)
            HEl = get_conn_helement(nI, nJ, ilutI, ilutJ, ic, excit, tParity, hel_unused)
        end if

        if (present(ncon)) ncon = ncon + 1

    end subroutine generate_connection_normal

    subroutine generate_connected_space_kpnt(original_space_size, original_space, &
                                             connected_space_size, connected_space, tSinglesOnlyOpt)

        use SymExcit2, only: gensymexcitit2par_worker

        ! This is the same as generate_connected_space, but using the old excitations
        ! generators because the new ones don't work with the tKPntSym option.

        ! This routine either counts or generates all the determinants connected to those in
        ! original_space. If connected_space is not present then they will only be counted,
        ! else they will be stored in connected_space. If tSinglesOnlyOpt is present and
        ! also .true. then actually this routine will only generate all single excitations
        ! (regardless of the system being studied), otherwise it will generate all connected
        ! determinants.

        use guga_bitRepOps, only: convert_ilut_toGUGA, convert_ilut_toNECI
        use guga_excitations, only: actHamiltonian
        use SystemData, only: tGUGA

        integer :: nexcit, j
        integer(n_int) :: ilutG(0:nifguga)
        integer(n_int), allocatable :: excitations(:, :)

        integer, intent(in) :: original_space_size
        integer(n_int), intent(in) :: original_space(0:, :)
        integer, intent(inout) :: connected_space_size
        integer(n_int), optional, intent(out) :: connected_space(0:, :)
        logical, intent(in), optional :: tSinglesOnlyOpt

        integer(n_int) :: ilutJ(0:NIfTot)
        integer :: nI(nel), nJ(nel)
        integer :: i, excit(2, 2), ex_flag
        integer :: nStore(6)
        integer, allocatable :: excit_gen(:)
        logical :: tStoreConnSpace, tSinglesOnly, tTempUseBrill, tAllExcitFound
        character(*), parameter :: this_routine = "generate_connected_space_kpnt"
        integer :: n_excits
        integer(n_int), allocatable :: temp_dets(:, :)

        if (present(connected_space)) then
            tStoreConnSpace = .true.
        else
            tStoreConnSpace = .false.
        end if

        tSinglesOnly = .false.
        if (present(tSinglesOnlyOpt)) then
            if (tSinglesOnlyOpt) tSinglesOnly = .true.
        end if

        connected_space_size = 0

        ! Over all the states in the original list:
        do i = 1, original_space_size

            call decode_bit_det(nI, original_space(0:NIfTot, i))

            if (t_k_space_hubbard) then

                ! for every loop we have to save the excitations per
                ! do we have to check if the list is unique?? i guess i do
#ifdef CMPLX_
                call stop_all(this_routine, "not implemented for complex")
#else
                call gen_all_excits_k_space_hubbard(nI, n_excits, temp_dets)
#endif

                if (tStoreConnSpace) then
                    connected_space(0:nifd, connected_space_size + 1:connected_space_size + n_excits) &
                        = temp_dets(0:nifd, :)
                end if
                connected_space_size = connected_space_size + n_excits

                ! GUGA changes:
                ! my actHamiltonian routine seems to satisfy the k-point symmetry
                ! so it should be straight forward to implemt it here too
            else if (tGUGA) then
                if (tSinglesOnly) call stop_all(this_routine, "don't use tSinglesOnly with GUGA")

                call convert_ilut_toGUGA(original_space(:, i), ilutG)

                call actHamiltonian(ilutG, CSF_Info_t(ilutG), excitations, nexcit)

                if (tStoreConnSpace) then
                    do j = 1, nexcit
                        call convert_ilut_toNECI(excitations(:, j), &
                                                 connected_space(:, connected_space_size + j))
                    end do
                end if

                connected_space_size = connected_space_size + nexcit

                deallocate(excitations)
                call LogMemDealloc(this_routine, tag_excitations)

            else

                call init_generate_connected_space(nI, ex_flag, tAllExcitFound, excit, excit_gen, nstore, tTempUseBrill)

                if (tSinglesOnly) ex_flag = 1

                do while (.true.)

                    call generate_connection_kpnt(nI, original_space(:, i), nJ, &
                                                  ilutJ, ex_flag, tAllExcitFound, nStore, excit_gen, &
                                                  ncon=connected_space_size)

                    if (tStoreConnSpace) then
                        connected_space(0:NIfD, connected_space_size) = ilutJ(0:NIfD)
                    end if

                end do

                call generate_connection_kpnt(nI, original_space(:, i), nJ, ilutJ, ex_flag, tAllExcitFound, &
                                              nStore, excit_gen, ncon=connected_space_size)

                if (tAllExcitFound) exit

                if (tStoreConnSpace) connected_space(0:NIfD, connected_space_size) = ilutJ(0:NIfD)

                deallocate(excit_gen)

            end if
        end do

    end subroutine generate_connected_space_kpnt

    subroutine generate_connection_kpnt(nI, ilutI, nJ, ilutJ, ex_flag, tAllExcitFound, nStore, &
                                        excit_gen, hel, ncon)

        use Determinants, only: get_helement
        use SymExcit2, only: gensymexcitit2par_worker

        integer :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel)
        integer(n_int), intent(out) :: ilutJ(0:NIfTot)
        integer, intent(inout) :: ex_flag
        logical, intent(inout) :: tAllExcitFound
        integer, intent(inout) :: nStore(6)
        integer, intent(inout) :: excit_gen(:)
        HElement_t(dp), optional, intent(out) :: hel
        integer, optional, intent(inout) :: ncon

        integer(n_int) :: ilut_tmp(0:NIfTot)
        integer :: ic

        ! Generate the next determinant.
        call GenSymExcitIt2Par_worker(nI, nel, G1, nBasis, .false., excit_gen, nJ, &
                                      ic, nStore, ex_flag, 1, nel)
        if (nJ(1) == 0) then
            tAllExcitFound = .true.
            return
        end if

        ! Encode nJ in ilut.
        call EncodeBitDet(nJ, ilutJ)

        ! If using HPHFs, and if this isn't the allowed state for this HPHF, find
        ! the the allowed state and continue with this as ilut.
        if (tHPHF .and. (.not. TestClosedShellDet(ilutJ))) then
            if (.not. IsAllowedHPHF(ilutJ(0:NIfD))) then
                ilut_tmp = ilutJ
                call FindExcitBitDetSym(ilut_tmp, ilutJ)
            end if
        end if

        if (present(HEl)) HEl = get_helement(nI, nJ, ic, ilutI, ilutJ)

        if (present(ncon)) ncon = ncon + 1

    end subroutine generate_connection_kpnt

end module enumerate_excitations
