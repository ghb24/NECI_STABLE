#include "macros.h"

! This module contains subroutines which generate all excitations from a given
! determinant, up to a certain excitation level. They output excitations one
! at a time, so that each time the subroutine is called the next excitation
! is generated. This is currently used by the semi-stochastic code to
! generate all states in the deterministic space.
module enumerate_excitations

    use bit_reps, only: NIfTot
    use constants
    use FciMCData, only: SpinInvBRR
    use Parallel_neci, only: MPISumAll
    use sort_mod, only: sort
    use SymData, only: nSymLabels
    use SymExcitDataMod
    use sym_general_mod
    use SystemData, only: nel, nBasis, G1, tFixLz, Arr, Brr

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
        integer :: ms_combination, ml_combination
        integer :: ms_num_combinations, ml_num_combinations
        integer :: total_ms, total_ml
        integer :: ms_1, ms_2, ml_1, ml_2
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
        integer :: orbI, ms_1_2, ml_1

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
                if (tFixLz) then
                    ml_1 = G1(orb1)%Ml
                else
                    ml_1 = 0
                end if
                ms_1_2 = 1 + iand(1, orb1)
                sym_ind1 = ClassCountInd(ms_1_2, G1(orb1)%Sym%S, ml_1)
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
        integer, pointer :: ms_combination_number, ml_combination_number
        integer, pointer :: ml_1, ml_2, total_ml, ml_num_combinations
        integer, pointer :: ml_combination
        integer :: min_ml, max_ml
        integer :: orbI, e1, e2, i1, i2, orb1a, orb2a, orbJ, lz_1, lz_2

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
        total_ml => gen_store%total_ml;  ml_num_combinations => gen_store%&
                                                       &ml_num_combinations
        ml_1 => gen_store%ml_1;          ml_2 => gen_store%ml_2
        ml_combination => gen_store%ml_combination

        ! Initialise the generator.
        if (ilut_ret(0) == -1) then
            sym1 = -1
            i = 1
            j = 0
            ms_combination = 1
            ml_combination = 1
        endif

        if (tFixLz) then
            min_ml = minval(G1(:)%Ml)
            max_ml = maxval(G1(:)%Ml)
        end if

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

                ! Find the number of different pairs of ms values that two electrons
                ! can have such that their total ms value is total_ms:
                call find_number_summation_possibilities(total_ms, -1, 1, &
                    ms_num_combinations, .true.)

                ! If applying Lz symmetry, find the number of pairs of ml values so
                ! that the total ml value is total_ml, as for ms above.
                if (tFixLz) then
                    total_ml = G1(orb1)%Ml + G1(orb2)%Ml
                    call find_number_summation_possibilities(total_ml, min_ml,&
                        max_ml, ml_num_combinations, .false.)
                else
                    ml_num_combinations = 1
                end if

                sym1 = 0
            end if

            do sym1 = sym1, nSymLabels - 1

                ! Generate the paired symmetry.
                sym2 = ieor(sym_prod, sym1)
                ! Don't overcount.
                if (sym2 < sym1) cycle

                do ms_combination = ms_combination, ms_num_combinations

                    do ml_combination = ml_combination, ml_num_combinations

                        ! If on the first loop of both this loop and the parent loop.
                        if (j == 0) then

                            ! Find the next pair of ms values that add up to total_ms.
                            call find_next_summation_pair(total_ms, -1, 1, ms_combination, &
                                .true., ms_1, ms_2)

                            ! Convert the ms values from +1 and -1 to 1 and 2 for alpha and
                            ! beta, respectively. This is required for ClassCoundInd.
                            ms_1 = (-ms_1+3)/2
                            ms_2 = (-ms_2+3)/2

                            ! If applying Lz symmetry, find the next pair of ml values that
                            ! add up to total_ml.
                            if (tFixLz) then
                                call find_next_summation_pair(total_ml, min_ml, max_ml, &
                                    ml_combination, .false., ml_1, ml_2)
                            else
                                ml_1 = 0
                                ml_2 = 0
                            end if

                            ! Finally, find the symmetry index for this specific (ms, sym, ml)
                            ! combination, for both electrons.
                            sym_ind1 = ClassCountInd(ms_1, sym1, ml_1)
                            sym_ind2 = ClassCountInd(ms_2, sym2, ml_2)
                            ind1 = SymLabelCounts2(1, sym_ind1)
                            ind2 = SymLabelCounts2(1, sym_ind2)
                            count1 = OrbClasscount(sym_ind1)
                            count2 = OrbClasscount(sym_ind2)
                            ! The total number of pairs of orbitals with these symmetries.
                            norb_sym = count1 * count2
                        end if

                        j = j + 1
                        ! Loop over all possble pairs of orbitals with the symmetries generated above.
                        do j = j, norb_sym

                            ! Direct mapping to orbitals.
                            i1 = mod(j-1, count1)
                            i2 = (j-1-i1)/count1
                            orbI = SymLabelList2(ind1 + i1)
                            orbJ = SymLabelList2(ind2 + i2)

                            ! If the two symmetries are the same (not including Ml or Ms)
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
                    ml_combination = 1

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

    subroutine find_number_summation_possibilities(target_integer, min_integer, max_integer, &
        num_pairs, ms_values)

        ! This routine finds the number of pairs of integers which sum to some given
        ! value, given that the two integers have a minimum and maximum size. If ms_values is
        ! true then the only the integers corresponding to actual ms values are cycled through.

        integer, intent(in) :: target_integer, min_integer, max_integer
        logical, intent(in) :: ms_values
        integer, intent(out) :: num_pairs
        integer :: i, j
        integer :: interval

        num_pairs = 0

        if (ms_values) then
            interval = 2
        else
            interval = 1
        end if

        ! Loop over all values for the first integer.
        do i = min_integer, max_integer, interval
            ! Loop over all values for the second integer.
            do j = min_integer, max_integer
                ! If the two integers add up to the target value.
                if ((i+j) == target_integer) num_pairs = num_pairs + 1
            end do
        end do

    end subroutine find_number_summation_possibilities

    subroutine find_next_summation_pair(target_integer, min_integer, max_integer, pair_number, &
        ms_values, number_1, number_2)

        ! This routine finds the next pair of numbers which add up to target_integer, given
        ! the maximum and minimum values of the induvidual numbers are min_integer and
        ! max_integer. pair_number specifies which generated pair should be returned.
        ! If ms_values is true then the only the integers corresponding to actual ms values
        ! are cycled through.

        integer, intent(in) :: target_integer, min_integer, max_integer, pair_number
        logical, intent(in) :: ms_values
        integer, intent(out) :: number_1, number_2
        integer :: i, j, counter
        integer :: interval = 1

        counter = 0

        if (ms_values) then
            interval = 2
        else
            interval = 1
        end if

        ! Loop over all values for the first integer.
        do i = min_integer, max_integer, interval
            ! Loop over all values for the second integer.
            do j = min_integer, max_integer
                ! If the two integers add up to the target value.
                if ((i+j) == target_integer) then
                    counter = counter+1
                    ! If on the pair specified by the user.
                    if (counter == pair_number) then
                        number_1 = i
                        number_2 = j
                        return
                    end if
                end if
            end do
        end do

    end subroutine find_next_summation_pair

    subroutine remove_high_energy_orbs(ilut_list, num_states, target_num_states, tParallel)

        integer, intent(inout) :: num_states
        integer(n_int), intent(inout) :: ilut_list(0:NIfTot, 1:num_states)
        integer, intent(in) :: target_num_states
        logical, intent(in) :: tParallel
        integer, allocatable, dimension(:) :: orbitals_rmvd
        integer :: i, j, orb, num_orbs_rmvd
        integer :: bit, elem
        integer :: states_rmvd_this_proc, counter
        integer(MPIArg) :: tot_num_states, states_rmvd_all_procs
        logical :: occupied

        states_rmvd_this_proc = 0

        if (tParallel) then
            call MPISumAll(int(num_states, MPIArg), tot_num_states)
        else
            tot_num_states = num_states
        end if

        if (tot_num_states <= target_num_states) return

        write(6,'(a32)') "Removing high energy orbitals..."

        ! Loop through all orbitals, from highest energy to lowest.
        do i = nBasis, 1, -1

            num_orbs_rmvd = nBasis-i+1

            orb = BRR(i)

            bit = mod(orb-1, bits_n_int)
            elem = (orb-1-bit)/bits_n_int

            ! Loop through all states and remove those states with orbital orb&
            ! occupied.
            do j = 1, num_states
                occupied = btest(ilut_list(elem, j), bit)
                if (occupied) then
                    ilut_list(:, j) = 0
                    states_rmvd_this_proc = states_rmvd_this_proc + 1
                end if
            end do

            if (tParallel) then
                call MPISumAll(int(states_rmvd_this_proc, MPIArg), states_rmvd_all_procs)
            else
                states_rmvd_all_procs = states_rmvd_this_proc
            end if

            ! If there are degenerate orbitals, then cycle to remove the
            ! degenrate orbitals too, before giving the program chance to quit.
            if (i > 1) then
                ! If the orbitals energies are the same:
                if ( Arr(i, 1) == Arr((i-1), 1) ) cycle
            end if

            if (tot_num_states-states_rmvd_all_procs <= target_num_states) exit

        end do

        ! Loop through the list and shuffle states down to fill in the gaps
        ! created above.
        counter = 0
        do i = 1, num_states
            ! If the state wasn't set to 0:
            if (.not. all(ilut_list(:,i) == 0)) then
                counter = counter + 1
                ilut_list(:, counter) = ilut_list(:, i)
            end if
        end do
        if (counter /= num_states-states_rmvd_this_proc) &
            call stop_all("remove_high_energy_orbs", &
                          "Wrong number of states found.")
        num_states = counter

        ! Finally, output information:
        write(6,'(a36)',advance='no') "The following orbitals were removed:"
        allocate(orbitals_rmvd(num_orbs_rmvd))
        orbitals_rmvd = BRR(nBasis-num_orbs_rmvd+1:nBasis)
        call sort(orbitals_rmvd)
        do i = 1, num_orbs_rmvd
            write(6,'(i5)',advance='no') orbitals_rmvd(i)
            if (i == num_orbs_rmvd) write(6,'()',advance='yes')
        end do
        deallocate(orbitals_rmvd)
        write(6,'(i6,1x,a29)') states_rmvd_all_procs, "states were removed in total."
        write(6,'(i6,1x,a17)') tot_num_states-states_rmvd_all_procs, "states were kept."
        call neci_flush(6)

    end subroutine remove_high_energy_orbs

end module enumerate_excitations
