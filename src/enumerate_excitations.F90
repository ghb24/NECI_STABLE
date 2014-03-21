#include "macros.h"

module enumerate_excitations

    use bit_rep_data, only: NIfD, NIfTot
    use bit_reps, only: decode_bit_det
    use constants
    use DetBitOps, only: IsAllowedHPHF, TestClosedShellDet, EncodeBitDet
    use FciMCData, only: SpinInvBRR
    use HPHFRandExcitMod, only: FindExcitBitDetSym
    use Parallel_neci, only: MPISumAll
    use sort_mod, only: sort
    use SymData, only: nSymLabels, SymTable, SymLabels, SymClasses
    use SymExcit3, only: GenExcitations3
    use SymExcitDataMod
    use sym_general_mod
    use SystemData, only: nel, nBasis, G1, tFixLz, Arr, Brr, tHPHF, tHub, &
                          tUEG, tKPntSym, tReal, tUseBrillouin

    implicit none

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
    
contains

    subroutine enumerate_spatial_excitations (ilutI, nI, ilut_ret, exflag, &
                                             gen_store)

        ! Generate all single and/or double spatial excitations (i.e. orbital
        ! configurations) of a given determinant. This is done so that, if an
        ! orbital is singly occupied, the electron is always placed in the beta
        ! orbital, as is the convention used for specifying CSFs. The spin
        ! eigenfunctions associated with an orbital configuration are *not*
        ! generated.

        ! In:  ilutI, nI - The determinant to excite
        !      exflag    - Set bits 0, 1 to indicate if singles/doubles should
        !                  be created
        ! IO:  gen_store - Stores the state of the generator
        ! Out: ilut_ret  - Returns the determinants produced. ilut_ret(0)
        !                  should be set to -1 to initialise, and will return
        !                  this once generation is complete.

        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer(n_int), intent(inout) :: ilut_ret(0:NIfTot)
        integer, intent(in) :: nI(nel), exflag
        logical :: bSingle, bDouble
        type(excit_store), intent(inout), target :: gen_store

        integer, pointer :: i, j, orb1, orb2, ind1, ind2, sym1, sym2
        integer, pointer :: count1, count2, sym_ind1, sym_ind2, norb_sym
        integer, pointer :: sym_prod, npairs
        integer :: orbI, orbJ, orb1a, orb2a, e1, e2, i1, i2

        ! Ensure we only generate the correct excitations
        bSingle = .false.
        bDouble = .false.
        if (btest(exflag, 0)) bSingle = .true.
        if (btest(exflag, 1)) bDouble = .true.

        ! Map the local variables onto the store
        i => gen_store%i;                j => gen_store%j
        orb1 => gen_store%orb1;          orb2 => gen_store%orb2
        ind1 => gen_store%ind1;          ind2 => gen_store%ind2
        sym1 => gen_store%sym1;          sym2 => gen_store%sym2
        count1 => gen_store%count1;      count2 => gen_store%count2
        sym_ind1 => gen_store%sym_ind1;  sym_ind2 => gen_store%sym_ind2
        norb_sym => gen_store%norb_sym
        sym_prod => gen_store%sym_prod
        npairs => gen_store%npairs

        ! Initialise the generator
        if (ilut_ret(0) == -1) then
            gen_store%gen_singles = .true.
            sym1 = -1
            i = 1
            j = 0
        endif

        ! Consider single excitations
        if (bSingle .and. gen_store%gen_singles) then
            ! Find the next possible single excitation. Loop over 
            ! electrons, and vacant orbitals. Interrupt loop when we
            ! find what we need.
            do i = i, nel

                if (j == 0) then
                    orb1 = nI(i)
                    sym_ind1 = ClassCountInd(2, G1(orb1)%Sym%S, 0)
                    orb2 = ab_pair(orb1)

                    ! We only want to excite one of doubly occupied pairs.
                    if (IsOcc(ilutI, orb2) .and. is_beta(orb1)) cycle

                    ind1 = SymLabelCounts2(1, sym_ind1)
                    norb_sym = OrbClasscount(sym_ind1)
                end if

                j = j + 1
                do j = j, norb_sym

                    orbI = SymLabelList2(ind1 + j - 1)

                    ! Cannot excite to the source of the excitation
                    if (is_in_pair(orb1, orbI)) cycle

                    ! Cannot excite to a doubly occupied orbital. If 
                    ! singly occupied, can only excite to the unoccupied
                    ! bit!
                    if (IsOcc(ilutI, orbI)) then
                        orbI = ab_pair(orbI)
                        if (IsOcc(ilutI, orbI)) cycle
                    endif

                    ! Now we can generate the determinant, and
                    ! interrupt the loop!!!!!
                    ilut_ret = ilutI
                    clr_orb(ilut_ret, orb1)
                    set_orb(ilut_ret, orbI)
                    return
                end do
                j = 0 ! Reset loop
            end do
            i = 1 ! Reset loop

            ! Finished generating singles
            gen_store%gen_singles = .false.
        else
            ! If we aren't generating singles, skip them.
            gen_store%gen_singles = .false.
        endif

        ! Consider double excitations
        if (bDouble) then

            npairs = nel * (nel - 1) / 2
            do i = i, npairs

                if (sym1 == -1) then
                    ! Pick electrons uniformly
                    e1 = ceiling((1.0_dp + sqrt(real(1 + 8*i,8))) / 2)
                    e2 = i - ((e1 - 1) * (e1 - 2)) / 2
                    orb1 = nI(e1);    orb1a = ab_pair(orb1)
                    orb2 = nI(e2);    orb2a = ab_pair(orb2)

                    ! If either of these is the beta electron from a doubly
                    ! occupied pair, and they are not from the same pair, 
                    ! cycle
                    if ( ((IsOcc(ilutI, orb1a) .and. is_beta(orb1)) .or.&
                          (IsOcc(ilutI, orb2a) .and. is_beta(orb2))) .and. &
                         .not. is_in_pair(orb1, orb2)) then
                        cycle
                    end if

                    sym_prod = ieor(G1(orb1)%Sym%S, G1(orb2)%Sym%S)

                    sym1 = 0
                end if

                do sym1 = sym1, nSymLabels - 1

                    if (j == 0) then
                        ! Generate the paired symmetry without double counting
                        sym2 = ieor(sym_prod, sym1)
                        if (sym2 < sym1) cycle

                        ! Symmetry indices
                        sym_ind1 = ClassCountInd(2, sym1, 0)
                        sym_ind2 = ClassCountInd(2, sym2, 0)
                        ind1 = SymLabelCounts2(1, sym_ind1)
                        ind2 = SymLabelCounts2(1, sym_ind2)
                        count1 = OrbClasscount(sym_ind1)
                        count2 = OrbClasscount(sym_ind2)
                        norb_sym = count1 * count2
                    end if

                    j = j + 1
                    do j = j, norb_sym

                        ! Direct mapping to orbitals
                        i1 = mod(j-1, count1)
                        i2 = (j-1-i1)/count1
                        orbI = SymLabelList2(ind1 + i1)
                        orbJ = SymLabelList2(ind2 + i2)

                        ! If the two symmetries are the same, only generate 
                        ! each pair one way around...
                        if (sym1 == sym2 .and. orbJ < orbI) cycle

                        ! We don't want to put two electrons into the same
                        ! spin orbital...
                        if (orbI == orbJ) orbJ = ab_pair(orbJ)

                        ! Exclude excitations to the source
                        if (is_in_pair(orbI, orb1) .or. &
                            is_in_pair(orbI, orb2)) cycle
                        if (is_in_pair(orbJ, orb1) .or. &
                            is_in_pair(orbJ, orb2)) cycle

                        ! Cannot excite to doubly occupied orbitals.
                        if (IsOcc(ilutI, orbI)) then
                            orbI = ab_pair(orbI)
                            if (IsOcc(ilutI, orbI) .or. &
                                is_in_pair(orbI, orbJ)) cycle
                        end if
                        if (IsOcc(ilutI, orbJ)) then
                            orbJ = ab_pair(orbJ)
                            if (IsOcc(ilutI, orbJ) .or. &
                                is_in_pair(orbI, orbJ)) cycle
                        end if

                        ! Now we can generate the determinant, and interrupt
                        ! the loop
                        ilut_ret = ilutI
                        clr_orb(ilut_ret, orb1)
                        clr_orb(ilut_ret, orb2)
                        set_orb(ilut_ret, orbI)
                        set_orb(ilut_ret, orbJ)
                        return
                    end do
                    j = 0 ! Break loop


                enddo
                sym1 = -1 ! Break loop

            end do

            ! And we are done!
            ilut_ret(0) = -1

        else
            ! If we aren't generating doubles, then we are done.
            ilut_ret(0) = -1
        endif

    end subroutine enumerate_spatial_excitations

    subroutine generate_connected_space(original_space_size, original_space, &
            connected_space_size, connected_space, tSinglesOnlyOpt)

        ! A wrapper function to call the correct generation routine.

        integer, intent(in) :: original_space_size
        integer(n_int), intent(in) :: original_space(0:NIfTot, original_space_size)
        integer, intent(inout) :: connected_space_size
        integer(n_int), optional, intent(out) :: connected_space(0:NIfTot, connected_space_size)
        logical, intent(in), optional :: tSinglesOnlyOpt

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

        ! This routine either counts or generates all the determinants connected to those in
        ! original_space. If connected_space is not present then they will only be counted,
        ! else they will be stored in connected_space. If tSinglesOnlyOpt is present and
        ! also .true. then actually this routine will only generate all single excitations
        ! (regardless of the system being studied), otherwise it will generate all connected
        ! determinants.

        integer, intent(in) :: original_space_size
        integer(n_int), intent(in) :: original_space(0:NIfTot, original_space_size)
        integer, intent(inout) :: connected_space_size
        integer(n_int), optional, intent(out) :: connected_space(0:NIfTot, connected_space_size)
        logical, intent(in), optional :: tSinglesOnlyOpt

        integer(n_int) :: ilut(0:NIfTot), ilut_tmp(0:NIfTot)
        integer :: nI(nel), nJ(nel)
        integer :: i, excit(2,2), ex_flag, ex_flag_temp
        logical :: tAllExcitFound, tParity, tStoreConnSpace, tSinglesOnly

        if (present(connected_space)) then
            tStoreConnSpace = .true.
        else
            tStoreConnSpace = .false.
        end if

        tSinglesOnly = .false.
        if (present(tSinglesOnlyOpt)) then
            if (tSinglesOnlyOpt) tSinglesOnly = .true.
        end if

        ! Which excitations levels should be considered?
        ! Singles only (1), doubles only (2) or singles and doubles (3).
        if (tSinglesOnly) then
            ex_flag = 1
        else
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
        end if

        connected_space_size = 0

        ! Over all the states in the original list:
        do i = 1, original_space_size

            call decode_bit_det(nI, original_space(:,i))

            tAllExcitFound = .false.
            excit = 0
            ex_flag_temp = ex_flag

            do while(.true.)

                ! Generate the next determinant.
                call GenExcitations3(nI, original_space(:,i), nJ, ex_flag_temp, excit, tParity, &
                                      tAllExcitFound, .false.)
                if (tAllExcitFound) exit

                ! Encode nJ in ilut.
                call EncodeBitDet(nJ, ilut)

                ! If using HPHFs, and if this isn't the allowed state for this HPHF, find
                ! the the allowed state and continue with this as ilut.
                if (tHPHF .and. (.not. TestClosedShellDet(ilut))) then
                    if (.not. IsAllowedHPHF(ilut(0:NIfD))) then
                        ilut_tmp = ilut
                        call FindExcitBitDetSym(ilut_tmp, ilut)
                    end if
                end if

                connected_space_size = connected_space_size + 1

                if (tStoreConnSpace) connected_space(0:NIfD, connected_space_size) = ilut(0:NIfD)

            end do

        end do

    end subroutine generate_connected_space_normal

    subroutine generate_connected_space_kpnt(original_space_size, original_space, &
            connected_space_size, connected_space, tSinglesOnlyOpt)

        ! This is the same as generate_connected_space, but using the old excitations
        ! generators because the new ones don't work with the tKPntSym option.

        ! This routine either counts or generates all the determinants connected to those in
        ! original_space. If connected_space is not present then they will only be counted,
        ! else they will be stored in connected_space. If tSinglesOnlyOpt is present and
        ! also .true. then actually this routine will only generate all single excitations
        ! (regardless of the system being studied), otherwise it will generate all connected
        ! determinants.

        integer, intent(in) :: original_space_size
        integer(n_int), intent(in) :: original_space(0:NIfTot, original_space_size)
        integer, intent(inout) :: connected_space_size
        integer(n_int), optional, intent(out) :: connected_space(0:NIfTot, connected_space_size)
        logical, intent(in), optional :: tSinglesOnlyOpt

        integer(n_int) :: ilut(0:NIfTot), ilut_tmp(0:NIfTot)
        integer :: nI(nel), nJ(nel)
        integer :: i, excit(2,2), ex_flag
        integer :: nStore(6), nExcitMemLen(1), iExcit, iMaxExcit, ierr
        integer, allocatable :: excit_gen(:)
        logical :: tStoreConnSpace, tSinglesOnly, tTempUseBrill
        character(*), parameter :: t_r = 'generate_connected_space_kpnt'

        if (present(connected_space)) then
            tStoreConnSpace = .true.
        else
            tStoreConnSpace = .false.
        end if

        tSinglesOnly = .false.
        if (present(tSinglesOnlyOpt)) then
            if (tSinglesOnlyOpt) tSinglesOnly = .true.
        end if

        ! Which excitations levels should be considered?
        ! Singles only (1), doubles only (2) or singles and doubles (3).
        if (tSinglesOnly) then
            ex_flag = 1
        else
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
        end if

        ! We have to ensure that brillouins theorem isn't on for the excitation
        ! generator.
        if (tUseBrillouin) then
            tTempUseBrill = .true.
            tUseBrillouin = .false.
        else
            tTempUseBrill = .false.
        end if

        connected_space_size = 0

        ! Over all the states in the original list:
        do i = 1, original_space_size

            call decode_bit_det(nI, original_space(:,i))

            ! This is all initialisation stuff for the excitation generator:

            iMaxExcit = 0
            nStore(1:6) = 0

            call GenSymExcitIt2(nI, nel, G1, nBasis, .true., nExcitMemLen, nJ, iMaxExcit, &
                                 nStore, ex_flag)

            allocate(excit_gen(nExcitMemLen(1)), stat=ierr)
            if (ierr .ne. 0) call Stop_All(t_r, "Problem allocating excitation generator")
            excit_gen = 0

            call GenSymExcitIt2(nI, nel, G1, nBasis, .true., excit_gen, nJ, iMaxExcit, &
                                 nStore, ex_flag)

            ! Initialisation finished.

            do while(.true.)

                ! Generate the next determinant.
                call GenSymExcitIt2(nI, nel, G1, nBasis, .false., excit_gen, nJ, iExcit, &
                                     nStore, ex_flag)
                if (nJ(1) == 0) exit

                ! Encode nJ in ilut.
                call EncodeBitDet(nJ, ilut)

                ! If using HPHFs, and if this isn't the allowed state for this HPHF, find
                ! the the allowed state and continue with this as ilut.
                if (tHPHF .and. (.not. TestClosedShellDet(ilut))) then
                    if (.not. IsAllowedHPHF(ilut(0:NIfD))) then
                        ilut_tmp = ilut
                        call FindExcitBitDetSym(ilut_tmp, ilut)
                    end if
                end if

                connected_space_size = connected_space_size + 1

                if (tStoreConnSpace) connected_space(0:NIfD, connected_space_size) = ilut(0:NIfD)

            end do

            deallocate(excit_gen)

        end do

    end subroutine generate_connected_space_kpnt

end module enumerate_excitations
