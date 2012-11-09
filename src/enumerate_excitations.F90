#include "macros.h"

! This module contains subroutines which generate all excitations from a given
! determinant, up to a certain excitation level. They output excitations one
! at a time, so that each time the subroutine is called the next excitation
! is generated. This is currently used by the semi-stochastic code to
! generate all states in the deterministic space.
module enumerate_excitations

    use bit_reps, only: NIfTot
    use constants
    use SymData, only: nSymLabels
    use SymExcitDataMod
    use sym_general_mod
    use SystemData, only: nel, G1

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
        integer :: ms_combination
        integer :: ms_num_combinations
        integer :: total_ms
        integer :: ms_1, ms_2
        logical :: gen_singles
    end type
    
contains

    subroutine enumerate_all_single_excitations(ilutI, nI, ilut_ret, gen_store)
  
        ! This subroutine generates all possible single excitations from a given determinant,
        ! within symmetry and spin restrictions. It will output a single determinant each time
        ! it is called. It uses a type to store the previous state of the variables in the
        ! subroutine so that the next time it is called it will generate the next determinant.

        ! In: ilutI, nI - The determinant to excite from.
        ! IO: gen_store - Stores the state of the generator.
        ! Out: ilut_ret - Returns the determinants produced. ilut_ret(0)
        ! should be set to -1 to initialise, and will return
        ! this once generation is complete.

        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer(n_int), intent(inout) :: ilut_ret(0:NIfTot)
        integer, intent(in) :: nI(nel)
        type(excit_store), intent(inout), target :: gen_store

        integer, pointer :: i, j, orb1, ind1, norb_sym, sym_ind1
        integer :: orbI, ms_1_2

        ! Map the local variables onto the store.
        i => gen_store%i;                j => gen_store%j
        orb1 => gen_store%orb1;          ind1 => gen_store%ind1;
        sym_ind1 => gen_store%sym_ind1;  norb_sym => gen_store%norb_sym

        ! Initialise the generator.
        if (ilut_ret(0) == -1) then
            i = 1
            j = 0
        endif

        ! Find the next possible single excitation. Loop over 
        ! electrons, and vacant orbitals. Interrupt loop when we
        ! find what we need.
        do i = i, nel

            if (j == 0) then
                orb1 = nI(i)
                ! This gets Ms such that beta is represented by 1 and
                ! alpha by 2, as required.
                ms_1_2 = 1 + iand(1, orb1)
                sym_ind1 = ClassCountInd(ms_1_2, G1(orb1)%Sym%S, 0)
                ind1 = SymLabelCounts2(1, sym_ind1)
                norb_sym = OrbClasscount(sym_ind1)
            end if

            j = j + 1
            do j = j, norb_sym

                orbI = SymLabelList2(ind1 + j - 1)

                ! Cannot excite to an occupied orbital
                if (IsOcc(ilutI, orbI)) cycle

                ! Generate the determinant and interrupt the loop.
                ilut_ret = ilutI
                clr_orb(ilut_ret, orb1)
                set_orb(ilut_ret, orbI)
                return
            end do
            j = 0 ! Reset loop
        end do

        ilut_ret(0) = -1

    end subroutine enumerate_all_single_excitations

    subroutine enumerate_all_double_excitations(ilutI, nI, ilut_ret, gen_store)

        ! This subroutine generates all possible double excitations from a given determinant,
        ! within symmetry and spin restrictions. It will output a single determinant each time
        ! it is called. It uses a type to store the previous state of the variables in the
        ! subroutine so that the next time it is called it will generate the next determinant.

        ! In: ilutI, nI - The determinant to excite from
        ! IO: gen_store - Stores the state of the generator.
        ! Out: ilut_ret - Returns the determinants produced. ilut_ret(0)
        ! should be set to -1 to initialise, and will return
        ! this once generation is complete.

        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer(n_int), intent(inout) :: ilut_ret(0:NIfTot)
        integer, intent(in) :: nI(nel)
        type(excit_store), intent(inout), target :: gen_store

        integer, pointer :: i, j, orb1, ind1, norb_sym, sym_ind1
        integer, pointer :: ms_combination, orb2, ind2, sym1, sym2
        integer, pointer :: count1, count2, npairs, sym_ind2, sym_prod
        integer, pointer :: ms_num_combinations, total_ms, ms_1, ms_2
        integer :: orbI, e1, e2, i1, i2, orb1a, orb2a, orbJ

        ! Map the local variables onto the store
        i => gen_store%i;                j => gen_store%j
        orb1 => gen_store%orb1;          orb2 => gen_store%orb2
        ind1 => gen_store%ind1;          ind2 => gen_store%ind2
        sym1 => gen_store%sym1;          sym2 => gen_store%sym2
        count1 => gen_store%count1;      count2 => gen_store%count2
        sym_ind1 => gen_store%sym_ind1;  sym_ind2 => gen_store%sym_ind2
        norb_sym => gen_store%norb_sym;  ms_combination => gen_store%ms_&
                                                              &combination
        sym_prod => gen_store%sym_prod;  npairs => gen_store%npairs
        total_ms => gen_store%total_ms;  ms_num_combinations => gen_store%&
                                                       &ms_num_combinations
        ms_1 => gen_store%ms_1;          ms_2 => gen_store%ms_2

        ! Initialise the generator.
        if (ilut_ret(0) == -1) then
            sym1 = -1
            i = 1
            j = 0
            ms_combination = 1
        endif

        npairs = nel * (nel - 1) / 2
        do i = i, npairs

            if (sym1 == -1) then
                ! Pick electrons uniformly
                e1 = ceiling((1.0 + sqrt(real(1 + 8*i))) / 2)
                e2 = i - ((e1 - 1) * (e1 - 2)) / 2
                orb1 = nI(e1);
                orb2 = nI(e2);

                sym_prod = ieor(G1(orb1)%Sym%S, G1(orb2)%Sym%S)
                total_ms = G1(orb1)%Ms + G1(orb2)%Ms
                ! If total_ms = -2 or +2 then both spins have to be beta
                ! or both alpha respectively, so there is only one combination
                ! to consider. If total_ms = 0 then the spins excited to can
                ! be (alpha, beta) or (beta, alpha). This formula maps -2, 0
                ! and +2 to the values 1, 2 and 1, as required.
                ms_num_combinations = (2-abs(total_ms))/2 + 1

                sym1 = 0
            end if

            do sym1 = sym1, nSymLabels - 1

                ! If total_ms = 0, loop over both possibilities
                do ms_combination = ms_combination, ms_num_combinations

                    if (j == 0) then
                        ! Set the Ms values (1=alpha, 2=beta) for both electrons.
                        if (total_ms == -2) then
                            ms_1 = 2
                            ms_2 = 2
                        else if (total_ms == 2) then
                            ms_1 = 1
                            ms_2 = 1
                        else
                            ms_1 = ms_combination ! 1 or 2.
                            ms_2 = 3 - ms_combination ! 2 or 1.
                        end if
                        ! Generate the paired symmetry.
                        sym2 = ieor(sym_prod, sym1)
                        ! Don't overcount.
                        if (sym2 < sym1) cycle

                        ! Symmetry indices, now including Ms values.
                        sym_ind1 = ClassCountInd(ms_1, sym1, 0)
                        sym_ind2 = ClassCountInd(ms_2, sym2, 0)
                        ind1 = SymLabelCounts2(1, sym_ind1)
                        ind2 = SymLabelCounts2(1, sym_ind2)
                        count1 = OrbClasscount(sym_ind1)
                        count2 = OrbClasscount(sym_ind2)
                        norb_sym = count1 * count2
                    end if

                    j = j + 1
                    do j = j, norb_sym

                        ! Direct mapping to orbitals.
                        i1 = mod(j-1, count1)
                        i2 = floor(real(j-1)/count1)
                        orbI = SymLabelList2(ind1 + i1)
                        orbJ = SymLabelList2(ind2 + i2)

                        ! If the two symmetries are the same (not including Ms),
                        ! only generate each pair one way around.
                        if (sym1 == sym2 .and. orbJ < orbI) cycle

                        ! We don't want to put two electrons into the same
                        ! spin orbital.
                        if (orbI == orbJ) cycle

                        ! Cannot excite to occupied orbitals.
                        if (IsOcc(ilutI, orbI)) cycle
                        if (IsOcc(ilutI, orbJ)) cycle

                        ! Now we can generate the determinant, and interrupt
                        ! the loop.
                        ilut_ret = ilutI
                        clr_orb(ilut_ret, orb1)
                        clr_orb(ilut_ret, orb2)
                        set_orb(ilut_ret, orbI)
                        set_orb(ilut_ret, orbJ)
                        return
                    end do
                    j = 0

                end do
                ms_combination = 1

            end do
            sym1 = -1

        end do

        ilut_ret(0) = -1

    end subroutine enumerate_all_double_excitations

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
                    e1 = ceiling((1.0 + sqrt(real(1 + 8*i))) / 2)
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
                        i2 = floor(real(j-1)/count1)
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

end module enumerate_excitations
