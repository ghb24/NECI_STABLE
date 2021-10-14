#include "macros.h"
! GUGA bit-representation operations:
! contains functions, concerning the representation of CSFs and
! associated information in NECI and operations thereon.
! some of these functions should be replaced by bit-ops defined in the
! macros.h header file, but for clarity and test purposes write them down
! explicetly too
module guga_bitRepOps

    use SystemData, only: nEl, Stot, nSpatOrbs, nbasis
    use guga_data, only: ExcitationInformation_t, excit_type, gen_type, &
                         rdm_ind_bitmask, pos_excit_lvl_bits, pos_excit_type_bits, &
                         n_excit_lvl_bits, n_excit_type_bits, n_excit_index_bits, &
                         excit_names
    use constants, only: dp, n_int, bits_n_int, bni_, bn2_, int_rdm, int64
    use DetBitOps, only: return_ms, count_set_bits, MaskAlpha, &
                         count_open_orbs, ilut_lt, ilut_gt, MaskAlpha, MaskBeta, &
                         CountBits, DetBitEQ
    use bit_rep_data, only: test_flag, flag_deltaB_single, IlutBits, &
                            flag_deltaB_double, flag_deltaB_sign, niftot, &
                            nIfGUGA, nIfd, BitRep_t, GugaBits
    use util_mod, only: binary_search, binary_search_custom, operator(.div.), &
                        near_zero, stop_all

    use sort_mod, only: sort

    use LoggingData, only: tRDMonfly

    use Fcimcdata, only: tFillingStochRDMonfly

    implicit none

    private
    public :: isDouble, csf_purify, get_preceeding_opposites, &
            isProperCSF_nI, isProperCSF_ilut, getDeltaB, setDeltaB, &
            encode_matrix_element, update_matrix_element, &
            extract_matrix_element, write_det_guga, &
            convert_ilut_toGUGA, convert_ilut_toNECI, convert_guga_to_ni, &
            write_guga_list, add_guga_lists, &
            findFirstSwitch, findLastSwitch, find_switches, &
            calcstepvector, &
            calcB_vector_int, calcB_vector_nI, calcB_vector_ilut, &
            count_open_orbs, count_open_orbs_ij, &
            count_beta_orbs_ij, count_alpha_orbs_ij, &
            calcOcc_vector_ilut, calcOcc_vector_int, &
            encodebitdet_guga, identify_excitation, &
            CSF_Info_t, current_csf_i, csf_ref, new_CSF_Info_t, fill_csf_i, &
            is_compatible, &
            calc_csf_i, extract_h_element, getexcitation_guga, &
            getspatialoccupation, getExcitationRangeMask, &
            contract_1_rdm_ind, contract_2_rdm_ind, extract_1_rdm_ind, &
            extract_2_rdm_ind, encode_rdm_ind, extract_rdm_ind, &
            encode_stochastic_rdm_x0, encode_stochastic_rdm_x1, &
            encode_stochastic_rdm_ind, encode_stochastic_rdm_info, &
            extract_stochastic_rdm_x0, extract_stochastic_rdm_x1, &
            extract_stochastic_rdm_ind, extract_stochastic_rdm_info, &
            init_guga_bitrep, transfer_stochastic_rdm_info, &
            extract_excit_lvl_rdm, extract_excit_type_rdm, &
            encode_excit_info, encode_excit_info_type, extract_excit_info_type, &
            encode_excit_info_indices, extract_excit_info_indices, &
            extract_excit_info, isProperCSF_flexible, find_guga_excit_lvl

    ! interfaces
    interface isProperCSF_ilut
        module procedure isProperCSF_b
        module procedure isProperCSF_sys
    end interface isProperCSF_ilut

    interface find_switches
        module procedure find_switches_ilut
        module procedure find_switches_stepvector
    end interface find_switches

    interface encode_matrix_element
        module procedure encode_matrix_element_real
#ifdef CMPLX_
        module procedure encode_matrix_element_cmplx
#endif
    end interface encode_matrix_element

    interface update_matrix_element
        module procedure update_matrix_element_real
#ifdef CMPLX_
        module procedure update_matrix_element_cmplx
#endif
    end interface update_matrix_element

    interface encode_excit_info
        module procedure encode_excit_info_scalar
        module procedure encode_excit_info_obj
        module procedure encode_excit_info_vec
    end interface encode_excit_info

    interface extract_excit_info_indices
        module procedure extract_excit_info_indices_vec
        module procedure extract_excit_info_indices_scalar
    end interface extract_excit_info_indices

    interface extract_excit_info
        module procedure extract_excit_info_scalar
        module procedure extract_excit_info_vector
        module procedure extract_excit_info_obj
    end interface extract_excit_info

    interface encode_excit_info_indices
        module procedure encode_excit_info_indices_vec
        module procedure encode_excit_info_indices_scalar
    end interface encode_excit_info_indices

    type :: CSF_Info_t
        integer, allocatable :: stepvector(:)
        integer, allocatable :: Occ_int(:), B_int(:)
        real(dp), allocatable :: Occ_real(:), B_real(:)

        real(dp), allocatable :: cum_list(:)
            !! also use a fake cum-list of the non-doubly occupied orbital to increase
            !! preformance in the picking of orbitals (a)
    end type

    interface CSF_Info_t
        module procedure construct_CSF_Info_t
    end interface

    type(CSF_Info_t) :: current_csf_i, csf_ref
        !! Information about the current CSF, similar to ilut and nI.
contains

    subroutine init_guga_bitrep(n_spatial_bits)
        ! set up a nIfGUGA variable to use a similar integer list to
        ! calculate excitations for a given GUGA CSF
        integer, intent(in) :: n_spatial_bits
        integer :: x0_pos, x1_pos, deltaB_pos, rdm_ind_pos, rdm_x0_pos, rdm_x1_pos

        ! Structure of a bit representation:
        ! the parenthesis is for the stochastic GUGA rdm implementation
        ! | 0-NIfD: Det | x0 | x1 | deltaB | (rdm_ind | rdm_x0 | rdm_x1)
        !
        ! -------
        ! (NIfD + 1) * 64-bits              Orbital rep.
        !  1         * 64-bits              x0 matrix element
        !  1         * 64-bits              x1 matrix element
        !  1         * 64-bits              deltaB value
        ! if we sample RDMs:
        !  1         * 64-bits              rdm_index (contains ex-level and type info!)
        !  1         * 64-bits              x0 coupling coeff for RDMs
        !  1         * 64-bits              x1-coupling coeff for RDMs

        x0_pos = n_spatial_bits + 1
        x1_pos = x0_pos + 1
        deltaB_pos = x1_pos + 1

        if (tRDMonfly) then
            rdm_ind_pos = deltaB_pos + 1
            rdm_x0_pos = rdm_ind_pos + 1
            rdm_x1_pos = rdm_x0_pos + 1

            nifguga = rdm_x1_pos
        else
            rdm_ind_pos = -1
            rdm_x0_pos = -1
            rdm_x1_pos = -1

            nIfGUGA = deltaB_pos
        end if

        ! and also use a global data structure for more overview
        GugaBits = BitRep_t( &
                   len_tot=nIfGUGA, &
                   len_orb=n_spatial_bits, &
                   ind_x0=x0_pos, &
                   ind_x1=x1_pos, &
                   ind_b=deltaB_pos, &
                   ind_rdm_ind=rdm_ind_pos, &
                   ind_rdm_x0=rdm_x0_pos, &
                   ind_rdm_x1=rdm_x1_pos)

    end subroutine init_guga_bitrep

    pure integer function find_guga_excit_lvl_to_doubles(ilutI, ilutJ)
        ! make an highly optimized excitation level finder for guga up
        ! to doubles (for efficiency reasons)
        ! maybe use this and the SD version to initialize a pointer
        ! at startup to use the correct one and to be sure to have a
        ! correct GUGA excit-lvl info at all necessary stages
        integer(n_int), intent(in) :: ilutI(0:GugaBits%len_tot), ilutJ(0:GugaBits%len_tot)

        unused_var(ilutI)
        unused_var(ilutJ)

        find_guga_excit_lvl_to_doubles = nel

        call stop_all("here", "todo")

    end function find_guga_excit_lvl_to_doubles


    pure integer function find_guga_excit_lvl(ilutI, ilutJ)
        ! general excit-level finder
        integer(n_int), intent(in) :: ilutI(0:GugaBits%len_tot), ilutJ(0:GugaBits%len_tot)

        type(ExcitationInformation_t) :: excitInfo

        excitInfo = identify_excitation(ilutI, ilutJ)

        find_guga_excit_lvl = excitInfo%excitLvl

    end function find_guga_excit_lvl

    ! finally with the, after all, usable trialz, popcnt, and leadz routines
    ! (except for the PGI and NAG compilers) i can write an efficient
    ! excitation identifier between two given CSFs
    ! for the big systems i realised that this is necessary, since
    ! applying the hamiltonian exactly to the reference derterminant is
    ! VERY time consuming!
    ! still this involves a lof of coding, since i do not have routines
    ! yet, which calculates the matrix element, if both CSFs are provided
    ! and probably also the whole excitation information can be provided

    pure function identify_excitation(ilutI, ilutJ) result(excitInfo)
        integer(n_int), intent(in) :: ilutI(0:nifd), ilutJ(0:nifd)
        type(ExcitationInformation_t) :: excitInfo

        integer(n_int) :: alpha_i(0:nifd), alpha_j(0:nifd), beta_i(0:nifd), &
                          beta_j(0:nifd), singles_i(0:nifd), singles_j(0:nifd), &
                          change_1(0:nifd), change_2(0:nifd), mask_singles(0:nifd), &
                          spin_change(0:nifd), overlap(0:nifd), &
                          mask_2(0:nifd), mask_3(0:nifd), mask_change_1(0:nifd), &
                          mask_change_2(0:nifd), mask_change_0(0:nifd)

        integer :: n_change_1, n_change_2, first_spin, last_spin, first_occ, &
                   last_occ, second_occ, third_occ, i, j, k, l, occ_double, &
                   ind(4), pos, ind_2(2), ind_3(2), res_orbs

        logical :: spin_change_flag

        ! i figure, that when i convert all the stepvectors like:
        ! 0 -> 0
        ! 1 -> 1
        ! 2 -> 1
        ! 3 -> 0
        ! which can be done by parts of the count open orbs routine
        ! i can identify the level of excitation(except for mixed full-starts
        ! full-stops..
        ! wait a minute ... what about 0 -> 3 type of excitations.. ??
        ! look at all the type of excitations:

        ! singles:
        ! di: 0123 -> 0110
        ! dj: 1212 -> 1111 -> 1001 -> correctly identifiable

        ! non-overlap:
        ! di: 0123 -> 0110
        ! dj: 1032 -> 1001 -> 1111 -> correct

        ! single overlap alike generators
        ! like singles and have to consider for identified singles also the
        ! matrix element influence of these.. but for the identification
        ! it is not a problem

        ! single overlap mixed generators:
        ! di: 1203 -> 1100
        ! dj: 0132 -> 0101 -> 1001 -> wrong!!
        ! have to identify the 0 -> 3 or 3 -> 0 switches indepentendly..
        ! also for the full-stop full-start type excitations with alike
        ! generators

        ! normal double:
        ! di: 0123 -> 0110
        ! dj: 1302 -> 1001 -> 1111 -> correct

        ! full-stop, or full-stop alike
        ! di: 0123 -> 0110
        ! dj: 1320 -> 1010 -> 1100 -> incorrect see above!

        ! full-stop or full-start mixed
        ! di: 1122 -> 1111
        ! dj: 1230 -> 1100 -> 0011 -> incorrect
        ! this is also a special case: it looks like a single, but i have
        ! to check if there is a change in the spin-coupling, below or
        ! above the identified single excitations

        ! full-start into full-stop alike:
        ! di: 0123 -> 0110
        ! dj: 3120 -> 0110 -> 0000 -> incorrect
        ! so i also have to check if it appears to be no excitation at all..
        ! puh ... thats gonna be tough..

        ! full-start into full-stop mixed:
        ! di: 1212 -> 1111
        ! dj: 1122 -> 1111 -> 0000 -> also looks like no excitation
        ! have to check for spin-coupling changes ..

        ! so:
        ! i have to look for occupation differences of +- 1
        ! i have to look for spin -> coupling changes
        ! i have to look for occupation differences of +- 2
        ! and i have to check the relation between the involved indices
        ! if it is a possible combination which leads to valid single
        ! or double exitations..

        ! the +-1 difference i figured out.

        ! what about spin-coupling changes?
        ! 123012
        ! 112120 i want:
        ! 010010 how?
        ! so i want the changes only in the singly occupied orbitals
        ! i could make a mask with the beginning stuff from above.. and then
        ! pick the xor only in these bits

        ! and i could also use the inversion of that mask to find the
        ! double occupation changes.. which i have to count twice anyway

        ! yeah this sounds not so bad.. but the painful part will be
        ! to find the relations of the hole and electron indices..
        ! but i have done that already kind of..

        ! AND i also have to write the matrix element calculation routine
        ! specifically for 2 given CSFs .. but atleast that is easier to do
        ! since i have all the necessary information at hand and dont need
        ! to do any branching and stuff.
        ! this could also make it easier to check if the excitation is
        ! compatible after all..

        ! and i have to set it up in such a way, that it can deal with
        ! multiple integers, for bigger number of orbitals, because thats
        ! exactly where it is necessary to use.
        ! i can use the stuff used in the paper by emmanuel to do that stuff
        ! i guess

        ! and since i need the created masks in all 3 of the checks i should
        ! not split up the calculation in multiple subroutine to be able
        ! to reuse them and be more efficient

        ! so the first part is to create a mask of the singly occupied
        ! orbitals and the inversion

        ! set the excit-level to some ridicolous high value for early returns
        excitInfo%excitLvl = nel
        excitInfo%valid = .false.
        ! some compilation problem again...
        ! first check if it is not the same ilut!
        if (DetBitEQ(ilutI, ilutJ)) then
            ! do diagonal element or return if i do not need the diagonals..
            excitInfo%excitLvl = 0
            return
        else

            ! i have to create the singles and non-singles mask and
            ! convert the 2 iluts to 1 -> 1, 2 -> 1, 0 -> 0, 3 -> 0
            alpha_i = iand(ilutI, MaskAlpha)
            beta_i = iand(ilutI, MaskBeta)
            alpha_i = ishft(alpha_i, -1)
            ! this is the ilut with 1 -> 1, 2 -> 1, 0 -> 0, 3 -> 0
            singles_i = ieor(alpha_i, beta_i)

            ! and do the same for J
            alpha_j = iand(ilutJ, MaskAlpha)
            beta_j = iand(ilutJ, MaskBeta)
            alpha_j = ishft(alpha_j, -1)

            singles_j = ieor(alpha_j, beta_j)

            ! with those 2 integers i can determine the +-1 occupation changes
            change_1 = ieor(singles_i, singles_j)
            ! can i tell something about the holes or electrons here

            ! i need the number and the indices maybe .. but i need the indices
            ! only if it is a valid single or double excitation, so do that
            ! later
            n_change_1 = sum(popcnt(change_1))
            ! remember a 2 would mean a valid single exitation here.. if no
            ! other stuff changed..

            ! and if i shift this to the right and xor with the original i get
            ! the singles mask
            ! the "correct" mask_singles:
            mask_singles = iand(singles_i, singles_j)
            mask_singles = ieor(mask_singles, ishft(mask_singles, +1))

            mask_change_1 = ieor(change_1, ishft(change_1, +1))
            ! and the doubles is the not of that

            ! actually for the double-change check i just need the change
            ! in the doubly occupied orbitals.. so where both alpha and beta
            ! orbitals are occupied

            ! the correct mask_doubles, atleast dn=+-2 or dn=0 & n=0,2
            mask_change_2 = iand(not(mask_change_1), not(mask_singles))
            mask_change_0 = iand(not(mask_change_1), mask_singles)

            ! also determine the +-2 changes
            ! should be enough to mask the positions of the 3 and 0 and then
            ! check if there are changes ..
            change_2 = ieor(iand(ilutI, mask_change_2), iand(ilutJ, mask_change_2))

            n_change_2 = sum(popcnt(change_2))
            ! can i say something about the holes and electrons here?
            ! i could check with the index of the change, if it is a 0 or a 3..

            ! so... what can i say about the validness of the excitation here ..
            ! i still have to determine spin coupling changes outside the
            ! index range of the excitation.. but it is probably better to do
            ! that later on and discard non-valid excitations at that point
            ! already..
            ! what is the maximum allowed number of changes?
            ! 4 i guess..
            ! i can have 4 singles changes -> non-overlap or proper double
            ! i can have 4 double changes -> full-start into fullstop alike
            ! i can have 2 singles and 2 doubles -> fullstart or fullstop alike
            ! and everything below ...
            ! and i can already include some decision making on what the
            ! type of excitations are which are possible..
            ! eg: n_change_2 = 4 & n_change_1 = 0 -> fullstart - fullstop alike
            ! if no other spin-coupling changes (inside and outside since x1 = 0)

            ! n_change_1 = 4 & n_change_2 = 0 -> non-overlap + normal double
            ! also if no other spin-change outside! the excitation range
            ! and matrix element depends on the spin-change inbetween
            ! which can exclude the non-overlap double!

            ! n_change_2 = 2 & n_change_1 = 0.. can that be? is this even a
            ! possibilit to detect? i dont think so.. this would mean a difference
            ! in electron number..

            ! n_change_1 = 2 & n_change_2 = 2 -> fullstart/fullstop
            ! + again checks outside the excitation range of spin-coupling changes

            ! n_change_1 = 2 & n_change_2 = 0 -> normal single or fullstart/fullstop
            ! mixed excitation. if only change above xor below excitation range
            ! not both! and i also have to consider all the possible singly
            ! occupied orbs which could lead to this excitation...

            ! n_change_1 = 0 & n_change 2 = 0 -> full-start into fullstop mixed
            ! this is the hard part.. here i have to consider all the possible
            ! index combination which could lead to this type of excitation

            ! but if n_change_1 + n_change_1 > 4 -> no possible excitation!

            ! hm i just realised in all the single changes <= 4 there is still
            ! a lot of other stuff that can happen..
            ! and i could also exit the routine if any of the other stuff
            ! does crazy shit, exept for the spin-coupling changes, which can
            ! happen a lot and still represent a valid excitation
            ! so atleast check the double changes too..

            if (n_change_1 + n_change_2 > 4) then
                ! no possible excitation
                return

            else
                ! could check spin-coupling here, since i need it in all cases
                ! below..
                spin_change = ieor(iand(ilutI, mask_singles), iand(ilutJ, mask_singles))

                ! need the first and last index to check the indices
                ! maybe i need to loop over the number of used integers to determine
                ! where the first and last spin change is in actual spatial
                ! orbitals..
                ! and i also need to consider that not all of the bits are used
                ! for the orbitals.. if not so many orbitals are involved..
                ! hm..
                ! in NECI the spin-orbitals get stored beginning from the left
                ! in the bit-representation.. so for leadz i have nothing to
                ! worry about the size of the basis
                ! what default value should i use to indicate no spin-coupling
                ! change?
                ! i only check if first spin is lower than anything.. so if i
                ! set it to nSpatOrbs + 1 it can never be true if no spin-coupling
                ! change is found
                first_spin = nSpatOrbs + 1
                ! i have to use bits_n_int not n_int!! to give me the number of
                ! orbitals in the integer!
                ! ok.. so now i know the convention in neci, it is actually
                ! stored in the "usual" way beginning from the right, but is
                ! just printed out differently in the end, which makes it a bit
                ! confusing..

                do i = 0, nifd
                    ! or should i exit
                    if (spin_change(i) == 0) cycle
                    ! so i have definetly something..
                    ! 1100 -> orbital 1 -> 1 + leadz=0 / 2
                    ! 0011 -> orbital 2 -> 1 + leadz=3 / 2
                    ! do i want it to be stored in spatial orbs already?
                    first_spin = 1 + ishft(bits_n_int * i + trailz(spin_change(i)), -1)

                    ! due to the if -> cycle statement i know that i have found
                    ! the FIRST spin change this way or none at all in the loop
                    ! so exit here
                    exit
                end do

                ! for the last_spin change i should default it to 0 so it is
                ! never true if no spin_change is ever found
                last_spin = 0
                ! maybe deal with the last integer of the bit representation
                ! specifically to identify the non-used part of the bits
                ! nah, with the modulo i do not need to do i specifically!
                ! yes do it specifically again! since so i can avoid doing the
                ! modulo all the time
                if (.not. spin_change(nifd) == 0) then
                    ! so i know it is in the last one, which i have to truncate
                    ! to the used bits

                    ! wtf... something wrong i think i need intermediate variable
                    ! due to some integer conversion stuff..
                    i = nBasis
                    j = leadz(spin_change(nifd)) - (bits_n_int - mod(nBasis, bits_n_int))

                    last_spin = ishft(i - j, -1)

                else
                    ! have to subtract the non-occupied orbs over the loop
                    ! and i already know the nidf-th had no change..
                    res_orbs = mod(nBasis, bits_n_int)

                    do i = nifd - 1, 0, -1
                        if (spin_change(i) == 0) then
                            ! have to increase res_orbs, but since here the hole
                            ! integer is used it is bits_n_int/2
                            res_orbs = res_orbs + bits_n_int
                            cycle
                        end if

                        ! just to be save from some integer conversion issues
                        j = nBasis
                        k = leadz(spin_change(i))

                        last_spin = ishft(j - res_orbs - k, -1)

                        ! then i can exit
                        exit
                    end do
                end if

                ! now i have the spin-changes in spatial orbitals already!!
                ! do that for the occupation indices also!

                ! so how do i deal with the case of multiple integers used for
                ! the spin-orbital storage??
                ! and also if not all bits are used in the integer representation!
                ! need to convert the obtained indices from each integer to
                ! the "actual" spin orbital(or maybe already to the spatial
                ! orbital)
                ! what if only 1 integer is used?
                ! and what is the convention in which order the orbitals are
                ! stored?? -> this decides when and where to use trailz or leadz

                ! hm.. i have to check for multiple integer use for more orbitals

                ! and determine the occupation indices also here .. nah, do not
                ! always need the same...

                ! do a select case, since more optimized
                select case (n_change_1)
                    ! what are the possible values? 0, 2, 4 i guess
                case (4)
                    ! this means n_change_2 = 0
                    ! and it is a non-overlap or "normal" double
                    ! still check spin-coupling outside excitation range

                    ! here i have to check if there is a spin-coupling change
                    ! outside the excitation range..
                    ! i know some changes happened!
                    ! i know here are 4 changes in change_1
                    j = 1
                    do i = 0, nifd
                        do while (change_1(i) /= 0)
                            pos = trailz(change_1(i))

                            ind(j) = 1 + ishft(bits_n_int * i + pos, -1)

                            ! i have to find out the convention in NECI how to
                            ! access the bits and how they are stored..
                            ! if they are really stored beginning from the left..
                            change_1(i) = ibclr(change_1(i), pos)

                            j = j + 1
                        end do
                    end do

                    first_occ = ind(1)
                    second_occ = ind(2)
                    third_occ = ind(3)
                    last_occ = ind(4)

                    ! i also have to deal with the edge-cases and if there are
                    ! no spin changes at all, what the compiler specific output
                    ! of leadz() and trailz() is for cases of no spin-changes
                    if (first_spin < first_occ .or. last_spin > last_occ) then
                        ! no excitation possible -> return

                        return
                    else
                        ! here the hard part of identifying the excitation
                        ! specifically comes in
                        ! how was it done in emmanuels paper to identify the
                        ! hole and electron?
                        ! i could use the isOne etc. macros with the indices
                        ! identified
                        ! determine the other indices..

                        ! i could also check if there is a spin-change in the
                        ! double overlap range, since if it is, there is no
                        ! possible non-overlap excitation
                        ! have to make a mask of ones in the overlap range
                        ! and then check with AND if there is any spin-change
                        ! in this range .. there are also some new fortran 2008
                        ! routines to create this masks

                        ! also the overlap has to be adapted to multiple
                        ! integer storage!
                        ! nah first i have to figure out the convention in neci
                        ! how the integer bits are accessed and outputted and
                        ! how the occupied orbitals are stored!
                        ! this messes up everything!
                        ! if the first index is in the 1. integer orbital
                        ! range
                        ! i need 2 indices to set both masks.. to which integer
                        ! the orbital index belongs to
                        ! mod(orb,bits_n_int) gives me the index in the integer
                        ! orb / bits_n_int, gives me the integer!
                        ! this does not work correctly yet!
                        ! i am not quite sure if this works as intended..
                        ! check!
                        ind_2 = [2 * second_occ / bits_n_int, mod(2 * second_occ, bits_n_int)]
                        ind_3 = [2 * (third_occ - 1) / bits_n_int, mod(2 * (third_occ - 1), bits_n_int)]
                        ! for the third index i have to put to 1 everyhing right
                        ! for the second, everything to the left
                        ! actually -1 has all bits set
                        mask_2(ind_2(1) + 1:nifd) = -1_n_int
                        mask_2(0:ind_2(1) - 1) = 0_n_int

                        mask_3(0:ind_3(1) - 1) = -1_n_int
                        mask_3(ind_3(1) + 1:nifd) = 0_n_int

                        ! so now in the mixed integer we have all 0
                        ! the maskl and maskr do not quite work as i suspected
                        ! and i have to work in spatial orbs! do not forget!
                        mask_2(ind_2(1)) = maskl(bits_n_int - ind_2(2), n_int)
                        mask_3(ind_3(1)) = maskr(ind_3(2), n_int)

                        overlap = iand(spin_change, iand(mask_2, mask_3))

                        if (sum(popcnt(overlap)) > 0) then
                            ! non-overlap NOT possible
                            spin_change_flag = .true.
                        else
                            spin_change_flag = .false.
                        end if

                        i = first_occ
                        j = second_occ
                        k = third_occ
                        l = last_occ

                        if (any(abs(calcB_vector_ilut(iLutI) - &
                            calcB_vector_ilut(ilutJ))>2)) return

                        ! puh.. this below MUST be optimized!!
                        ! but how?
                        ! how must i put in the order parameter.. to get the
                        ! correct sign in the 2-body integral contributions..
                        ! TODO i have discrepancies in the assignment of this
                        ! order signs.. in the nosym_ implementation i take care
                        ! of them and depending on the ordering of the generators
                        ! i assign +1 or -1, but in the sym_ implementation i
                        ! seem to totally neglect it..
                        ! did i update it in such a way that the sign-assignment
                        ! is not necessary anymore.. then i have to get rid of
                        ! it in the nosym_ approach too i guess, since i calculate
                        ! with only +1 signs afterwards, or did i jsut fuck it up
                        ! in the sym_ approach by totally neglecting them..
                        ! either way i think on of them is wrong..
                        ! welp, i definetly use it to get the matrix elements
                        ! at semi-starts and semi-stops, and also in the
                        ! 2-body contribution.. or did i find a way to
                        ! store the i,j,k,l indices in the sym_ approach in such
                        ! a way that this sign is alway +1?
                        ! but why are there no errors in the test-cases???
                        ! something is messed up with that... damn..
                        ! figure that out! that could lead to totally wrong matrix
                        ! elements.. atleast it should.. but why again does it
                        ! not show up in the tests??
                        ! for now stick to the Shavitt paper convention..
                        if (isZero(ilutI, i)) then
                            if (isZero(ilutI, j)) then
                                ! _R(i) -> _RR(j) -> ^RR(k) -> ^R(l)
                                excitInfo = assign_excitInfo_values_exact( &
                                            excit_type%double_raising, &
                                            gen_type%R, gen_type%R, gen_type%R, gen_type%R, gen_type%R, &
                                            j, l, i, k, i, j, k, l, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                            else if (isThree(ilutI, j)) then
                                ! have to check where the electron goes
                                if (isZero(ilutI, k)) then
                                    ! _R(i) -> _LR(j) -> ^LR(k) -> ^R(l)
                                    excitInfo = assign_excitInfo_values_exact( &
                                                excit_type%double_R_to_L_to_R, &
                                                gen_type%R, gen_type%L, gen_type%R, gen_type%R, gen_type%R, &
                                                i, l, k, j, i, j, k, l, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                                else if (isThree(ilutI, k)) then
                                    ! _R(i) -> _LR(j) -> ^RL(k) -> ^L(l)
                                    excitInfo = assign_excitInfo_values_exact( &
                                                excit_type%double_R_to_L, &
                                                gen_type%R, gen_type%L, gen_type%R, gen_type%R, gen_type%L, &
                                                i, k, l, j, i, j, k, l, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                                else
                                    ! n(k) = 1
                                    if (isZero(ilutJ, k)) then
                                        ! _R(i) -> _LR(j) -> ^RL(k) -> ^L(l)
                                        excitInfo = assign_excitInfo_values_exact( &
                                                    excit_type%double_R_to_L, &
                                                    gen_type%R, gen_type%L, gen_type%R, gen_type%R, gen_type%L, &
                                                    i, k, l, j, i, j, k, l, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                                    else
                                        ! _R(i) -> _LR(j) -> ^LR(k) -> ^R(l)
                                        excitInfo = assign_excitInfo_values_exact( &
                                                    excit_type%double_R_to_L_to_R, &
                                                    gen_type%R, gen_type%L, gen_type%R, gen_type%R, gen_type%R, &
                                                    i, l, k, j, i, j, k, l, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                                    end if
                                end if
                            else
                                ! n(j) = 1
                                if (isZero(ilutJ, j)) then
                                    ! _R(i) -> _LR(j) ...
                                    if (isZero(ilutI, k)) then
                                        ! _R(i) -> _LR(j) -> ^LR(k) -> ^R(l)
                                        excitInfo = assign_excitInfo_values_exact( &
                                                    excit_type%double_R_to_L_to_R, &
                                                    gen_type%R, gen_type%L, gen_type%R, gen_type%R, gen_type%R, &
                                                    i, l, k, j, i, j, k, l, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                                    else if (isThree(ilutI, k)) then
                                        ! _R(i) -> _LR(j) -> ^RL(k) -> ^L(l)
                                        excitInfo = assign_excitInfo_values_exact( &
                                                    excit_type%double_R_to_L, &
                                                    gen_type%R, gen_type%L, gen_type%R, gen_type%R, gen_type%L, &
                                                    i, k, l, j, i, j, k, l, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                                    else
                                        ! n(k) = 1
                                        if (isZero(ilutJ, k)) then
                                            ! _R(i) -> _LR(j) -> ^RL(k) -> ^L(l)
                                            excitInfo = assign_excitInfo_values_exact( &
                                                        excit_type%double_R_to_L, &
                                                        gen_type%R, gen_type%L, gen_type%R, gen_type%R, gen_type%L, &
                                                        i, k, l, j, i, j, k, l, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                                        else
                                            ! _R(i) -> _LR(j) -> ^LR(k) -> ^R(l)
                                            excitInfo = assign_excitInfo_values_exact( &
                                                        excit_type%double_R_to_L_to_R, &
                                                        gen_type%R, gen_type%L, gen_type%R, gen_type%R, gen_type%R, &
                                                        i, l, k, j, i, j, k, l, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                                        end if
                                    end if
                                else
                                    ! _R(i) -> _RR(j) -> ^RR(k) -> ^R(l)
                                    excitInfo = assign_excitInfo_values_exact( &
                                                excit_type%double_raising, &
                                                gen_type%R, gen_type%R, gen_type%R, gen_type%R, gen_type%R, &
                                                j, l, i, k, i, j, k, l, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)
                                end if
                            end if
                        else if (isThree(ilutI, i)) then
                            ! _L(i) -> ..
                            if (isZero(ilutI, j)) then
                                ! _L(i) -> _RL(j) -> ...
                                if (isZero(ilutI, k)) then
                                    ! _L(i) -> _RL(j) -> ^LR(k) -> ^R(l)
                                    excitInfo = assign_excitInfo_values_exact( &
                                                excit_type%double_L_to_R, &
                                                gen_type%R, gen_type%L, &
                                                gen_type%L, gen_type%L, gen_type%R, &
                                                j, l, k, i, i, j, k, l, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                                else if (isThree(ilutI, k)) then
                                    ! _L(i) -> _RL(j) -> ^RL(k) -> ^L(l)
                                    excitInfo = assign_excitInfo_values_exact( &
                                                excit_type%double_L_to_R_to_L, &
                                                gen_type%R, gen_type%L, &
                                                gen_type%L, gen_type%L, gen_type%L, &
                                                j, k, l, i, i, j, k, l, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                                else
                                    ! n(k) = 1
                                    if (isZero(ilutJ, k)) then
                                        ! _L(i) -> _RL(j) -> ^RL(k) -> ^L(l)
                                        excitInfo = assign_excitInfo_values_exact( &
                                                    excit_type%double_L_to_R_to_L, &
                                                    gen_type%R, gen_type%L, &
                                                    gen_type%L, gen_type%L, gen_type%L, &
                                                    j, k, l, i, i, j, k, l, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                                    else
                                        ! _L(i) -> _RL(j) -> ^LR(k) -> ^R(l)
                                        excitInfo = assign_excitInfo_values_exact( &
                                                    excit_type%double_L_to_R, &
                                                    gen_type%R, gen_type%L, &
                                                    gen_type%L, gen_type%L, gen_type%R, &
                                                    j, l, k, i, i, j, k, l, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                                    end if
                                end if
                            else if (isThree(ilutI, j)) then
                                ! _L(i) -> _LL(j) -> ^LL(k) -> ^L(l)
                                excitInfo = assign_excitInfo_values_exact( &
                                            excit_type%double_lowering, &
                                            gen_type%L, gen_type%L, &
                                            gen_type%L, gen_type%L, gen_type%L, &
                                            k, i, l, j, i, j, k, l, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                            else
                                ! n(j) = 1
                                if (isZero(ilutJ, j)) then
                                    ! _L(i) -> _LL(j) -> ^LL(k) -> ^L(l)
                                    excitInfo = assign_excitInfo_values_exact( &
                                                excit_type%double_lowering, &
                                                gen_type%L, gen_type%L, &
                                                gen_type%L, gen_type%L, gen_type%L, &
                                                k, i, l, j, i, j, k, l, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                                else
                                    ! _L(i) -> _RL(j) -> ...
                                    if (isZero(ilutI, k)) then
                                        ! _L(i) -> _RL(j) -> ^LR(k) -> ^R(l)
                                        excitInfo = assign_excitInfo_values_exact( &
                                                    excit_type%double_L_to_R, &
                                                    gen_type%R, gen_type%L, &
                                                    gen_type%L, gen_type%L, gen_type%R, &
                                                    j, l, k, i, i, j, k, l, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                                    else if (isThree(ilutI, k)) then
                                        ! _L(i) -> _RL(j) -> ^RL(k) -> ^L(l)
                                        excitInfo = assign_excitInfo_values_exact( &
                                                    excit_type%double_L_to_R_to_L, &
                                                    gen_type%R, gen_type%L, &
                                                    gen_type%L, gen_type%L, gen_type%L, &
                                                    j, k, l, i, i, j, k, l, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                                    else
                                        ! n(k) = 1
                                        if (isZero(ilutJ, k)) then
                                            ! _L(i) -> _RL(j) -> ^RL(k) -> ^L(l)
                                            excitInfo = assign_excitInfo_values_exact( &
                                                        excit_type%double_L_to_R_to_L, &
                                                        gen_type%R, gen_type%L, &
                                                        gen_type%L, gen_type%L, gen_type%L, &
                                                        j, k, l, i, i, j, k, l, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                                        else
                                            ! _L(i) -> _RL(j) -> ^LR(k) -> ^R(l)
                                            excitInfo = assign_excitInfo_values_exact( &
                                                        excit_type%double_L_to_R, &
                                                        gen_type%R, gen_type%L, &
                                                        gen_type%L, gen_type%L, gen_type%R, &
                                                        j, l, k, i, i, j, k, l, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                                        end if
                                    end if
                                end if
                            end if
                        else
                            ! n(i) = 1
                            if (isZero(ilutJ, i)) then
                                ! _L(i) -> ...
                                if (isZero(ilutI, j)) then
                                    ! _L(i) -> _RL(j) -> ...
                                    if (isZero(ilutI, k)) then
                                        ! _L(i) -> _RL(j) -> ^LR(k) -> ^R(l)
                                        excitInfo = assign_excitInfo_values_exact( &
                                                    excit_type%double_L_to_R, &
                                                    gen_type%R, gen_type%L, &
                                                    gen_type%L, gen_type%L, gen_type%R, &
                                                    j, l, k, i, i, j, k, l, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                                    else if (isThree(ilutI, k)) then
                                        ! _L(i) -> _RL(j) -> ^RL(k) -> ^L(l)
                                        excitInfo = assign_excitInfo_values_exact( &
                                                    excit_type%double_L_to_R_to_L, &
                                                    gen_type%R, gen_type%L, &
                                                    gen_type%L, gen_type%L, gen_type%L, &
                                                    j, k, l, i, i, j, k, l, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                                    else
                                        ! n(k) = 1
                                        if (isZero(ilutJ, k)) then
                                            ! _L(i) -> _RL(j) -> ^RL(k) -> ^L(l)
                                            excitInfo = assign_excitInfo_values_exact( &
                                                        excit_type%double_L_to_R_to_L, &
                                                        gen_type%R, gen_type%L, &
                                                        gen_type%L, gen_type%L, gen_type%L, &
                                                        j, k, l, i, i, j, k, l, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                                        else
                                            ! _L(i) -> _RL(j) -> ^LR(k) -> ^R(l)
                                            excitInfo = assign_excitInfo_values_exact( &
                                                        excit_type%double_L_to_R, &
                                                        gen_type%R, gen_type%L, &
                                                        gen_type%L, gen_type%L, gen_type%R, &
                                                        j, l, k, i, i, j, k, l, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                                        end if
                                    end if
                                else if (isThree(ilutI, j)) then
                                    ! _L(i) -> _LL(j) -> ^LL(k) -> ^L(l)
                                    excitInfo = assign_excitInfo_values_exact( &
                                                excit_type%double_lowering, &
                                                gen_type%L, gen_type%L, &
                                                gen_type%L, gen_type%L, gen_type%L, &
                                                k, i, l, j, i, j, k, l, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                                else
                                    ! n(j) = 1
                                    if (isZero(ilutJ, j)) then
                                        ! _L(i) -> _LL(j) -> ^LL(k) -> ^L(l)
                                        excitInfo = assign_excitInfo_values_exact( &
                                                    excit_type%double_lowering, &
                                                    gen_type%L, gen_type%L, &
                                                    gen_type%L, gen_type%L, gen_type%L, &
                                                    k, i, l, j, i, j, k, l, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                                    else
                                        ! _L(i) -> _RL(j) -> ...
                                        if (isZero(ilutI, k)) then
                                            ! _L(i) -> _RL(j) -> ^LR(k) -> ^R(l)
                                            excitInfo = assign_excitInfo_values_exact( &
                                                        excit_type%double_L_to_R, &
                                                        gen_type%R, gen_type%L, &
                                                        gen_type%L, gen_type%L, gen_type%R, &
                                                        j, l, k, i, i, j, k, l, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                                        else if (isThree(ilutI, k)) then
                                            ! _L(i) -> _RL(j) -> ^RL(k) -> ^L(l)
                                            excitInfo = assign_excitInfo_values_exact( &
                                                        excit_type%double_L_to_R_to_L, &
                                                        gen_type%R, gen_type%L, &
                                                        gen_type%L, gen_type%L, gen_type%L, &
                                                        j, k, l, i, i, j, k, l, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                                        else
                                            ! n(k) = 1
                                            if (isZero(ilutJ, k)) then
                                                ! _L(i) -> _RL(j) -> ^RL(k) -> ^L(l)
                                                excitInfo = assign_excitInfo_values_exact( &
                                                            excit_type%double_L_to_R_to_L, &
                                                            gen_type%R, gen_type%L, &
                                                            gen_type%L, gen_type%L, gen_type%L, &
                                                            j, k, l, i, i, j, k, l, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                                            else
                                                ! _L(i) -> _RL(j) -> ^LR(k) -> ^R(l)
                                                excitInfo = assign_excitInfo_values_exact( &
                                                            excit_type%double_L_to_R, &
                                                            gen_type%R, gen_type%L, &
                                                            gen_type%L, gen_type%L, gen_type%R, &
                                                            j, l, k, i, i, j, k, l, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                                            end if
                                        end if
                                    end if
                                end if
                            else
                                ! _R(i) -> ...
                                if (isZero(ilutI, j)) then
                                    ! _R(i) -> _RR(j) -> ^RR(k) -> ^R(l)
                                    excitInfo = assign_excitInfo_values_exact( &
                                                excit_type%double_raising, &
                                                gen_type%R, gen_type%R, &
                                                gen_type%R, gen_type%R, gen_type%R, &
                                                j, l, i, k, i, j, k, l, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                                else if (isThree(ilutI, j)) then
                                    ! have to check where the electron goes
                                    if (isZero(ilutI, k)) then
                                        ! _R(i) -> _LR(j) -> ^LR(k) -> ^R(l)
                                        excitInfo = assign_excitInfo_values_exact( &
                                                    excit_type%double_R_to_L_to_R, &
                                                    gen_type%R, gen_type%L, &
                                                    gen_type%R, gen_type%R, gen_type%R, &
                                                    i, l, k, j, i, j, k, l, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                                    else if (isThree(ilutI, k)) then
                                        ! _R(i) -> _LR(j) -> ^RL(k) -> ^L(l)
                                        excitInfo = assign_excitInfo_values_exact( &
                                                    excit_type%double_R_to_L, &
                                                    gen_type%R, gen_type%L, &
                                                    gen_type%R, gen_type%R, gen_type%L, &
                                                    i, k, l, j, i, j, k, l, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                                    else
                                        ! n(k) = 1
                                        if (isZero(ilutJ, k)) then
                                            ! _R(i) -> _LR(j) -> ^RL(k) -> ^L(l)
                                            excitInfo = assign_excitInfo_values_exact( &
                                                        excit_type%double_R_to_L, &
                                                        gen_type%R, gen_type%L, &
                                                        gen_type%R, gen_type%R, gen_type%L, &
                                                        i, k, l, j, i, j, k, l, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                                        else
                                            ! _R(i) -> _LR(j) -> ^LR(k) -> ^R(l)
                                            excitInfo = assign_excitInfo_values_exact( &
                                                        excit_type%double_R_to_L_to_R, &
                                                        gen_type%R, gen_type%L, &
                                                        gen_type%R, gen_type%R, gen_type%R, &
                                                        i, l, k, j, i, j, k, l, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                                        end if
                                    end if
                                else
                                    ! n(j) = 1
                                    if (isZero(ilutJ, j)) then
                                        ! _R(i) -> _LR(j) ...
                                        if (isZero(ilutI, k)) then
                                            ! _R(i) -> _LR(j) -> ^LR(k) -> ^R(l)
                                            excitInfo = assign_excitInfo_values_exact( &
                                                        excit_type%double_R_to_L_to_R, &
                                                        gen_type%R, gen_type%L, &
                                                        gen_type%R, gen_type%R, gen_type%R, &
                                                        i, l, k, j, i, j, k, l, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                                        else if (isThree(ilutI, k)) then
                                            ! _R(i) -> _LR(j) -> ^RL(k) -> ^L(l)
                                            excitInfo = assign_excitInfo_values_exact( &
                                                        excit_type%double_R_to_L, &
                                                        gen_type%R, gen_type%L, &
                                                        gen_type%R, gen_type%R, gen_type%L, &
                                                        i, k, l, j, i, j, k, l, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                                        else
                                            ! n(k) = 1
                                            if (isZero(ilutJ, k)) then
                                                ! _R(i) -> _LR(j) -> ^RL(k) -> ^L(l)
                                                excitInfo = assign_excitInfo_values_exact( &
                                                            excit_type%double_R_to_L, &
                                                            gen_type%R, gen_type%L, &
                                                            gen_type%R, gen_type%R, gen_type%L, &
                                                            i, k, l, j, i, j, k, l, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                                            else
                                                ! _R(i) -> _LR(j) -> ^LR(k) -> ^R(l)
                                                excitInfo = assign_excitInfo_values_exact( &
                                                            excit_type%double_R_to_L_to_R, &
                                                            gen_type%R, gen_type%L, &
                                                            gen_type%R, gen_type%R, gen_type%R, &
                                                            i, l, k, j, i, j, k, l, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                                            end if
                                        end if
                                    else
                                        ! _R(i) -> _RR(j) -> ^RR(k) -> ^R(l)
                                        excitInfo = assign_excitInfo_values_exact( &
                                                    excit_type%double_raising, &
                                                    gen_type%R, gen_type%R, &
                                                    gen_type%R, gen_type%R, gen_type%R, &
                                                    j, l, i, k, i, j, k, l, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)
                                    end if
                                end if
                            end if
                        end if
                    end if

                case (2)
                    ! n_change_2 could still be 0 or 2

                    select case (n_change_2)

                    case (2)
                        ! this means a a fullstart/fullstop still need to
                        ! check spin-coupling
                        ! can also be a single overlap with mixed generators!

                        j = 1
                        do i = 0, nifd
                            do while (change_1(i) /= 0)
                                pos = trailz(change_1(i))

                                ind(j) = 1 + ishft(bits_n_int * i + pos, -1)

                                change_1(i) = ibclr(change_1(i), pos)

                                j = j + 1
                            end do
                        end do
                        first_occ = ind(1)
                        last_occ = ind(2)

                        ! there is only one +-2 change!
                        do i = 0, nifd
                            if (change_2(i) == 0) cycle

                            occ_double = 1 + ishft(bits_n_int * i + trailz(change_2(i)), -1)

                            exit
                        end do

                        if (first_spin < min(first_occ, occ_double) .or. &
                            last_spin > max(last_occ, occ_double)) then
                            ! no excitation possible

                            return

                        else
                            ! identify the excitation specifically
                            ! there is no spin-coupling change allowed in the
                            ! double excitation region, since x1 matrix element
                            ! is 0 in this case anyway, remember that!!

                            ! have to check if its a single overlap, to check
                            ! for spin coupling changes within DE range
                            if (occ_double < first_occ) then
                                ! it is a fullstart alike -> check spin-coupling
                                ! i already know the first spin change is not
                                ! totally out of the excitation range!
                                if (first_spin < first_occ) then
                                    ! thats all i need to check, since first_spin
                                    ! is definetly lower than last_spin and i just
                                    ! need to check if there is any spin-coupling
                                    ! change in the DE region

                                    ! no excitation possible!

                                    return

                                else
                                    ! it is a full-start alike!
                                    ! have to still check if there is some
                                    ! spin-coupling change in the overlap region
                                    ! if there is -> no valid excitation!
                                    ind_2 = [2 * occ_double / bits_n_int, mod(2 * occ_double, bits_n_int)]
                                    ind_3 = [2 * (first_occ - 1) / bits_n_int, mod(2 * (first_occ - 1), bits_n_int)]
                                    ! for the third index i have to put to 1 everyhing right
                                    ! for the second, everything to the left
                                    mask_2(ind_2(1) + 1:nifd) = -1_n_int
                                    mask_2(0:ind_2(1) - 1) = 0_n_int

                                    mask_3(0:ind_3(1) - 1) = -1_n_int
                                    mask_3(ind_3(1) + 1:nifd) = 0_n_int

                                    ! so now in the mixed integer we have all 0
                                    mask_2(ind_2(1)) = maskl(bits_n_int - ind_2(2), n_int)
                                    mask_3(ind_3(1)) = maskr(ind_3(2), n_int)

                                    overlap = iand(spin_change, iand(mask_2, mask_3))

                                    if (sum(popcnt(overlap)) > 0) then
                                        ! no excitation possible!

                                        return

                                    else
                                        ! valid excitation
                                        i = occ_double
                                        j = first_occ
                                        k = last_occ

                                        if (any(abs(calcB_vector_ilut(iLutI) - &
                                            calcB_vector_ilut(ilutJ))>2)) return

                                        if (isZero(ilutI, i)) then
                                            ! _RR_(i) -> ^RR(j) -> ^R(k)
                                            excitInfo = assign_excitInfo_values_exact( &
                                                        excit_type%fullstart_raising, &
                                                        gen_type%R, gen_type%R, &
                                                        gen_type%R, gen_type%R, gen_type%R, &
                                                        i, j, i, k, i, i, j, k, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                                        else
                                            ! _LL_(i) -> ^LL(j) -> ^L(k)
                                            excitInfo = assign_excitInfo_values_exact( &
                                                        excit_type%fullstart_lowering, &
                                                        gen_type%L, gen_type%L, &
                                                        gen_type%L, gen_type%L, gen_type%L, &
                                                        k, i, j, i, i, i, j, k, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                                        end if
                                    end if
                                end if

                            else if (occ_double > last_occ) then
                                ! fullstop alike -> check if last spin-coupling
                                ! change is in the DE range
                                if (last_spin > last_occ) then
                                    ! no excitation possible

                                    return

                                else
                                    ! it is a valid fullstop alike excitation
                                    ! still have to check spin-coupling change
                                    ind_2 = [2 * last_occ / bits_n_int, mod(2 * last_occ, bits_n_int)]
                                    ind_3 = [2 * (occ_double - 1) / bits_n_int, mod(2 * (occ_double - 1), bits_n_int)]
                                    ! for the third index i have to put to 1 everyhing right
                                    ! for the second, everything to the left
                                    mask_2(ind_2(1) + 1:nifd) = -1_n_int
                                    mask_2(0:ind_2(1) - 1) = 0_n_int

                                    mask_3(0:ind_3(1) - 1) = -1_n_int
                                    mask_3(ind_3(1) + 1:nifd) = 0_n_int

                                    ! so now in the mixed integer we have all 0
                                    mask_2(ind_2(1)) = maskl(bits_n_int - ind_2(2), n_int)
                                    mask_3(ind_3(1)) = maskr(ind_3(2), n_int)

                                    overlap = iand(spin_change, iand(mask_2, mask_3))

                                    if (sum(popcnt(overlap)) > 0) then
                                        ! no excitation possible

                                        return

                                    else
                                        ! valid!
                                        i = first_occ
                                        j = last_occ
                                        k = occ_double

                                        if (any(abs(calcB_vector_ilut(iLutI) - &
                                            calcB_vector_ilut(ilutJ))>2)) return

                                        if (isZero(ilutI, k)) then
                                            ! _L(i) -> _LL(j) -> ^LL^(k)
                                            excitInfo = assign_excitInfo_values_exact( &
                                                        excit_type%fullstop_lowering, &
                                                        gen_type%L, gen_type%L, &
                                                        gen_type%L, gen_type%L, gen_type%L, &
                                                        k, i, k, j, i, j, k, k, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                                        else
                                            ! _R(i) -> _RR(j) -> ^RR^(k)
                                            excitInfo = assign_excitInfo_values_exact( &
                                                        excit_type%fullstop_raising, &
                                                        gen_type%R, gen_type%R, &
                                                        gen_type%R, gen_type%R, gen_type%R, &
                                                        i, k, j, k, i, j, k, k, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                                        end if
                                    end if
                                end if

                            else
                                ! the double change is between the single
                                ! changes -> so it is a single overlap mixed
                                ! no need to check the spin-coupling
                                i = first_occ
                                j = occ_double
                                k = last_occ

                                if (any(abs(calcB_vector_ilut(iLutI) - &
                                    calcB_vector_ilut(ilutJ))>2)) return

                                if (isZero(ilutI, j)) then
                                    ! _L(i) > ^LR_(j) -> ^R(k)
                                    excitInfo = assign_excitInfo_values_exact( &
                                                excit_type%single_overlap_L_to_R, &
                                                gen_type%R, gen_type%L, &
                                                gen_type%L, gen_type%L, gen_type%R, &
                                                j, k, j, i, i, j, j, k, 0, 2, 1.0_dp, 1.0_dp, 1, spin_change_flag)

                                else
                                    ! _R(i) -> ^RL_(j) -> ^L(k)
                                    excitInfo = assign_excitInfo_values_exact( &
                                                excit_type%single_overlap_R_to_L, &
                                                gen_type%R, gen_type%L, &
                                                gen_type%R, gen_type%R, gen_type%L, &
                                                i, j, k, j, i, j, j, k, 0, 2, 1.0_dp, 1.0_dp, 1, spin_change_flag)

                                end if
                            end if
                        end if

                    case (0)
                        ! "normal" single -> check spin-coupling

                        ! in this case a spin change either below or above is
                        ! possible, but no both..
                        j = 1
                        do i = 0, nifd
                            do while (change_1(i) /= 0)
                                pos = trailz(change_1(i))

                                ind(j) = 1 + ishft(bits_n_int * i + pos, -1)

                                ! i have to find out the convention in NECI how to
                                ! access the bits and how they are stored..
                                ! if they are really stored beginning from the left..
                                change_1(i) = ibclr(change_1(i), pos)

                                j = j + 1
                            end do
                        end do

                        first_occ = ind(1)
                        last_occ = ind(2)

                        if (first_spin < first_occ .and. last_spin > last_occ) then
                            ! no excitation possible

                            return

                        else if (first_spin < first_occ) then
                            ! full-start mixed
                            i = first_spin
                            j = first_occ
                            k = last_occ

                            if (any(abs(calcB_vector_ilut(iLutI) - &
                                calcB_vector_ilut(ilutJ))>2)) return

                            if (isZero(ilutI, j)) then
                                ! _RL_(i) -> ^LR(j) -> ^R(k)
                                excitInfo = assign_excitInfo_values_exact( &
                                            excit_type%fullstart_L_to_R, &
                                            gen_type%R, gen_type%L, &
                                            gen_type%R, gen_type%R, gen_type%R, &
                                            i, k, j, i, i, i, j, k, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                            else if (isThree(ilutI, j)) then
                                ! _RL_(i) -> ^RL(j) -> ^L(k)
                                excitInfo = assign_excitInfo_values_exact( &
                                            excit_type%fullstart_R_to_L, &
                                            gen_type%R, gen_type%L, &
                                            gen_type%R, gen_type%R, gen_type%L, &
                                            i, j, k, i, i, i, j, k, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                            else
                                if (isZero(ilutJ, j)) then
                                    ! _RL_(i) -> ^RL(j) -> ^L(k)
                                    excitInfo = assign_excitInfo_values_exact( &
                                                excit_type%fullstart_R_to_L, &
                                                gen_type%R, gen_type%L, &
                                                gen_type%R, gen_type%R, gen_type%L, &
                                                i, j, k, i, i, i, j, k, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                                else
                                    ! _RL_(i) -> ^LR(j) -> ^R(k)
                                    excitInfo = assign_excitInfo_values_exact( &
                                                excit_type%fullstart_L_to_R, &
                                                gen_type%R, gen_type%L, &
                                                gen_type%R, gen_type%R, gen_type%R, &
                                                i, k, j, i, i, i, j, k, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                                end if
                            end if

                        else if (last_spin > last_occ) then
                            ! full-stop mixed

                            i = first_occ
                            j = last_occ
                            k = last_spin

                            if (any(abs(calcB_vector_ilut(iLutI) - &
                                calcB_vector_ilut(ilutJ))>2)) return

                            if (isZero(ilutI, i)) then
                                ! _R(i) -> _LR(j) -> ^RL^(k)
                                excitInfo = assign_excitInfo_values_exact( &
                                            excit_type%fullstop_R_to_L, &
                                            gen_type%R, gen_type%L, &
                                            gen_type%R, gen_type%R, gen_type%R, &
                                            i, k, k, j, i, j, k, k, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                            else if (isThree(ilutI, i)) then
                                ! _L(i) -> _RL(j) -> ^RL^(k)
                                excitInfo = assign_excitInfo_values_exact( &
                                            excit_type%fullstop_L_to_R, &
                                            gen_type%R, gen_type%L, &
                                            gen_type%L, gen_type%L, gen_type%R, &
                                            j, k, k, i, i, j, k, k, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                            else
                                if (isZero(ilutJ, i)) then
                                    ! _L(i) -> _RL(j) -> ^RL^(k)
                                    excitInfo = assign_excitInfo_values_exact( &
                                                excit_type%fullstop_L_to_R, &
                                                gen_type%R, gen_type%L, &
                                                gen_type%L, gen_type%L, gen_type%R, &
                                                j, k, k, i, i, j, k, k, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                                else
                                    ! _R(i) -> _LR(j) -> ^RL^(k)
                                    excitInfo = assign_excitInfo_values_exact( &
                                                excit_type%fullstop_R_to_L, &
                                                gen_type%R, gen_type%L, &
                                                gen_type%R, gen_type%R, gen_type%R, &
                                                i, k, k, j, i, j, k, k, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                                end if
                            end if

                        else
                            ! regular single excitation
                            ! this is the "easiest" case.. start with this ..
                            ! 0123 -> 0110
                            ! 1212 -> 1111 -> 1001

                            ! 0123 -> 0110
                            ! 0312 -> 0011 -> 0101

                            ! 0123 -> 0110
                            ! 0303 -> 0000 -> 0110
                            ! -> the original stepvector number on the identified
                            ! indices is not enough to determine if it is a
                            ! hole or electron, if they are singly occupied ...

                            ! first i have to change the spin-orbital indices
                            ! to spatial orbitals: use beta orbs by definition
                            i = first_occ
                            j = last_occ

                            if (any(abs(calcB_vector_ilut(iLutI) - &
                                calcB_vector_ilut(ilutJ))>1)) return

                            ! but remember i need to calculate ALL the possible
                            ! double excitation influences which can lead
                            ! to this single excitation!
                            if (isZero(ilutI, i)) then
                                ! then i know first_occ is a hole and the second
                                ! a electron..
                                ! -> so it is a raising generator
                                ! _R(i) -> ^R(j)
                                excitInfo = assign_excitInfo_values_single_ex(gen_type%R, i, j, i, j)

                            else if (isThree(ilutI, i)) then
                                ! i know (i) is a electron ->
                                ! _L(i) -> ^L(j)
                                excitInfo = assign_excitInfo_values_single_ex(gen_type%L, j, i, i, j)

                            else
                                ! i know n(i) = 1
                                ! so i only need to check what happened in ilutJ
                                if (isZero(ilutJ, i)) then
                                    ! (i) was an electron
                                    ! _L(i) -> ^L(j)
                                    excitInfo = assign_excitInfo_values_single_ex(gen_type%L, j, i, i, j)

                                else
                                    ! (i) was a hole
                                    ! _R(i) -> ^R(j)
                                    excitInfo = assign_excitInfo_values_single_ex(gen_type%R, i, j, i, j)

                                end if
                            end if
                        end if
                    end select

                case (0)
                    ! n_change_2 can be 0 or 4, since 2 is not possible to
                    ! conserve electron number
                    select case (n_change_2)

                    case (4)
                        ! fullstart -> fullstop alike, check spin-coupling
                        ! due to x1 element being zero in this case there is
                        ! no spin-coupling change at all allowed..
                        if (sum(popcnt(spin_change)) > 0) then
                            ! no excitation possible

                            return

                        else

                            j = 1
                            do i = 0, nifd
                                do while (change_2(i) /= 0)
                                    pos = trailz(change_2(i))

                                    ind(j) = 1 + ishft(bits_n_int * i + pos, -1)

                                    ! i have to find out the convention in NECI how to
                                    ! access the bits and how they are stored..
                                    ! if they are really stored beginning from the left..
                                    change_2(i) = ibclr(change_2(i), pos)
                                    ! i need to clear the 2 set bits in the change 2
                                    change_2(i) = ibclr(change_2(i), pos + 1)

                                    j = j + 1
                                end do
                            end do

                            first_occ = ind(1)
                            last_occ = ind(2)
                            ! here excitation identification should be not so
                            ! hard .. but eg. here there should be now spin-change
                            ! also within in the excitation range, since the
                            ! x1-matrix element branch is 0...
                            ! also above for the fullstart/fullstop alikes..
                            ! so actually there is no spin-coupling change allowed
                            ! in this case..

                            i = first_occ
                            j = last_occ

                            if (any(abs(calcB_vector_ilut(iLutI) - &
                                calcB_vector_ilut(ilutJ))>2)) return

                            if (isZero(ilutI, i)) then
                                ! _RR_(i) -> ^RR^(j)
                                excitInfo = assign_excitInfo_values_exact( &
                                            excit_type%fullstart_stop_alike, &
                                            gen_type%R, gen_type%R, &
                                            gen_type%R, gen_type%R, gen_type%R, &
                                            i, j, i, j, i, i, j, j, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                            else
                                ! _LL_(i) -> ^LL^(j)
                                excitInfo = assign_excitInfo_values_exact( &
                                            excit_type%fullstart_stop_alike, &
                                            gen_type%L, gen_type%L, &
                                            gen_type%L, gen_type%L, gen_type%L, &
                                            j, i, j, i, i, i, j, j, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                            end if
                        end if

                    case (0)
                        ! full-start -> fullstop mixed or diagonal matrix element!
                        ! i could exclude diagonal from the beginning, if
                        ! no change at all between die iluts
                        ! have excluded diagonals above so only check for the
                        ! first and last spin-coupling changes

                        ! here there are a bunch of excitations possible..
                        ! maybe check for the first and last spin change and
                        ! then output a flag to calculate all the possible
                        ! excitations

                        ! check atleast the first and last spin-coupling changes..
                        i = first_spin
                        j = last_spin

                        if (any(abs(calcB_vector_ilut(iLutI) - &
                            calcB_vector_ilut(ilutJ))>2)) return

                        ! _RL_(i) -> ^RL^(j)
                        excitInfo = assign_excitInfo_values_exact( &
                                    excit_type%fullstart_stop_mixed, &
                                    gen_type%R, gen_type%L, gen_type%R, gen_type%R, gen_type%R, &
                                    i, j, j, i, i, i, j, j, 0, 2, 1.0_dp, 1.0_dp, 2, spin_change_flag)

                    end select
                end select
            end if
        end if

    end function identify_excitation

    pure function assign_excitInfo_values_exact(typ, gen1, gen2, currentGen, firstGen, &
                lastGen, i, j, k, l, fullStart, secondStart, firstEnd, fullEnd, &
                weight, excitLvl, order, order1, overlap, spin_change) &
                result(excitInfo)
        ! version of the excitation information filler for the exact
        ! matrix element calculation between 2 given CSFs
        integer, intent(in) :: typ, gen1, gen2, currentGen, firstGen, lastGen, &
                               i, j, k, l, fullStart, secondStart, firstEnd, &
                               fullEnd, weight, excitLvl
        integer, intent(in), optional :: overlap
        real(dp), intent(in) :: order, order1
        logical, intent(in), optional :: spin_change
        type(ExcitationInformation_t) :: excitInfo

        ! todo: asserts!
        excitInfo%typ = typ
        excitInfo%gen1 = gen1
        excitInfo%gen2 = gen2
        excitInfo%currentGen = currentGen
        excitInfo%firstGen = firstGen
        excitInfo%lastGen = lastGen
        excitInfo%i = i
        excitInfo%j = j
        excitInfo%k = k
        excitInfo%l = l
        excitInfo%fullStart = fullStart
        excitInfo%secondStart = secondStart
        excitInfo%firstEnd = firstEnd
        excitInfo%fullEnd = fullEnd
        excitInfo%weight = weight
        excitInfo%excitLvl = excitLvl
        excitInfo%order = order
        excitInfo%order1 = order1
        if (present(spin_change)) then
            excitInfo%spin_change = spin_change
        else
            excitInfo%spin_change = .false.
        end if
        if (present(overlap)) then
            excitInfo%overlap = overlap
        else
            excitInfo%overlap = 2
        end if

        excitInfo%valid = .true.

    end function assign_excitInfo_values_exact

    pure function assign_excitInfo_values_single_ex(gen, i, j, fullStart, fullEnd, typ) &
        result(excitInfo)
        integer, intent(in) :: gen, i, j, fullStart, fullEnd
        integer, intent(in), optional :: typ
        type(ExcitationInformation_t) :: excitInfo

        ! set default values for single excitations: which cause errors if
        ! called in incorrect places
        if (present(typ)) then
            excitInfo%typ = typ
        else
            excitInfo%typ = excit_type%single
        end if

        if (i == j) then
            excitInfo%excitLvl = 0
            excitInfo%weight = i
        else
            excitInfo%excitLvl = 1
            excitInfo%weight = 0
        end if
        excitInfo%k = 0
        excitInfo%l = 0
        excitInfo%secondStart = 0
        excitInfo%firstEnd = 0
        excitInfo%gen2 = -2
        excitInfo%order = 0.0_dp
        excitInfo%order1 = 0.0_dp
        excitInfo%overlap = 0

        ! then set proper values
        excitInfo%i = i
        excitInfo%j = j
        excitInfo%gen1 = gen
        excitInfo%fullStart = fullStart
        excitInfo%fullEnd = fullEnd
        excitInfo%currentGen = gen
        excitInfo%firstGen = gen
        excitInfo%lastGen = gen

        excitInfo%valid = .true.

    end function assign_excitInfo_values_single_ex

    subroutine getExcitation_guga(nI, nJ, ex)
        ! routine to determine excitation in guga basis
        ! for now to it very naively and unellegant by converting to ilut
        ! format and calculating occupation vectors
        integer, intent(in) :: nI(nEl), nJ(nEl)
        integer, intent(out) :: ex(2, 2)

        integer(n_int) :: ilutI(0:niftot), ilutJ(0:niftot)
        integer :: first, last, cnt_e, cnt_h, occ_diff(nSpatOrbs), i

        call EncodeBitDet_guga(nI, ilutI)
        call EncodeBitDet_guga(nJ, ilutJ)

        occ_diff = calcOcc_vector_int(ilutI(0:nifd)) - calcOcc_vector_int(ilutJ(0:nifd))

        select case (sum(abs(occ_diff)))

        case (0)
            ! there is a different definition of RDMs in the guga case or?
            ! because for _RL_ -> ^RL^ excitations or nI = nJ i end up here
            ! but a lot of orbital index combinations can lead to
            ! this type of excitation.. where should i assign it to..?
            ! read the R.Shephard paper and think of a new calculation of
            ! the RDM calculation in the GUGA case..

            ! ... for now, but in the indices of the first and last switch..

            first = findFirstSwitch(ilutI, ilutJ, 1, nSpatOrbs)

            last = findLastSwitch(ilutI, ilutJ, first, nSpatOrbs)

            ex(1, 1) = 2 * first
            ex(1, 2) = 2 * last - 1

            ex(2, 1) = 2 * first - 1
            ex(2, 2) = 2 * last

        case (2)

            ! this is a "normal" double excitation
            ! find the electron in nI which gets excited

            ex(1, 2) = 0
            ex(2, 2) = 0

            do i = 1, nSpatOrbs
                if (occ_diff(i) == 1) ex(1, 1) = 2 * i
                if (occ_diff(i) == -1) ex(2, 1) = 2 * i
            end do

        case (4)

            ! keep count of the already found electrons and holes
            cnt_e = 1
            cnt_h = 1

            do i = 1, nSpatOrbs

                select case (occ_diff(i))

                    ! this choice of default spin-orbitals below
                    ! makes certain two_rdm samplings as default alos..
                    ! not sure if this choice alone is valid..
                    ! also have to ensure i get the "spins" right so the
                    ! rest of the NECI RDM routines can handle that..
                case (2)
                    ! two eletrons get excited from orb i:
                    ex(1, 1) = 2 * i - 1
                    ex(1, 2) = 2 * i

                case (1)
                    ! one electron gets excited from i
                    ! at the first encountered electron cnt_e = 1
                    ! -> so this below gives me an alpha electron!
                    ! -> at the second it will be an beta to ensure
                    ! i get "correct" spins..
                    ex(1, cnt_e) = 2 * i - cnt_e + 1

                    cnt_e = cnt_e + 1

                case (-1)
                    ! one hole found
                    ex(2, cnt_h) = 2 * i - cnt_h + 1

                    cnt_h = cnt_h + 1

                case (-2)
                    ! two electron get excited to orb i
                    ex(2, 1) = 2 * i - 1
                    ex(2, 2) = 2 * i

                end select

            end do

        end select

    end subroutine getExcitation_guga

    subroutine find_switches_ilut(ilut, ind, lower, upper)
        ! for single excitations this checks for available switches around an
        ! already chosen index
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        integer, intent(in) :: ind
        integer, intent(out) :: lower, upper

        integer :: i
        ! set defaults if no such switches are available
        lower = 1
        upper = nSpatOrbs

        ! i think i could improve that with the new trailz, leadz routines..

        if (isOne(ilut, ind)) then
            do i = ind - 1, 2, -1
                if (isTwo(ilut, i)) then
                    lower = i
                    exit
                end if
            end do
            do i = ind + 1, nSpatOrbs - 1
                if (isTwo(ilut, i)) then
                    upper = i
                    exit
                end if
            end do
        else if (isTwo(ilut, ind)) then
            do i = ind - 1, 2, -1
                if (isOne(ilut, i)) then
                    lower = i
                    exit
                end if
            end do
            do i = ind + 1, nSpatOrbs - 1
                if (isOne(ilut, i)) then
                    upper = i
                    exit
                end if
            end do
        end if

    end subroutine find_switches_ilut

    subroutine find_switches_stepvector(csf_i, ind, lower, upper)
        ! same as above but using the already calculated stepvector
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: ind
        integer, intent(out) :: lower, upper
        character(*), parameter :: this_routine = "find_switches_stepvector"

        integer :: switch, i
        ! set defaults
        lower = 1
        upper = nSpatOrbs

        if (csf_i%stepvector(ind) == 1) then
            switch = 2

        else if (csf_i%stepvector(ind) == 2) then
            switch = 1

        else
            ! wrong input!
            call stop_all(this_routine, "wrong input! stepvalue /= {1,2}!")
        end if

        do i = ind - 1, 2, -1
            if (csf_i%stepvector(i) == switch) then
                lower = i
                exit
            end if
        end do
        do i = ind + 1, nSpatOrbs - 1
            if (csf_i%stepvector(i) == switch) then
                upper = i
                exit
            end if
        end do

    end subroutine find_switches_stepvector

    pure function findFirstSwitch(iI, iJ, start, semi) result(orb)
        ! write a scratch implementation to find the first change in
        ! stepvector for two given CSFs. do it inefficiently for now
        ! improve later on
        integer(n_int), intent(in) :: iI(0:GugaBits%len_tot), iJ(0:GugaBits%len_tot)
        integer, intent(in) :: start, semi
        integer :: orb, a, b

        integer :: i
        ! with the fortran 2008 intrinsic funcitons it would be easy...
        ! for now just do a loop over double overlap region and compare
        ! stepvalues

        ! i could also use the quantity here.. or?
        ! if it is always called for the current looked at ilut..
        ! i guess it does..

        ! implement this with the new f2008 routines..
        ! i need to find the first spin-change between start and semi-1

        if (start >= semi) then
            orb = -1
            return
        end if
        ! make the spin_change bit-rep

        ! todo check if this change worked!
        ! ok... i have two different goals here..
        ! before i wanted to check for any switches.. now i only want
        ! spin-changes.. to i ever need anything else then spin-changes?
        ! ok i really only need spin-changes.. so change the testsuite
        orb = -1
        do i = start, semi - 1
            a = getStepvalue(iI, i)
            b = getStepvalue(iJ, i)
            if (a /= b) then
                orb = i
                return
            end if
        end do

    end function findFirstSwitch

    pure function findLastSwitch(ilutI, ilutJ, semi, ende) result(orb)
        ! function to find last switch in a mixed fullstop excitation
        integer(n_int), intent(in) :: ilutI(0:GugaBits%len_tot), ilutJ(0:GugaBits%len_tot)
        integer, intent(in) :: ende, semi
        integer :: orb

        integer :: a, b, iOrb


        ! set it to impossible value, so contribution does not get
        ! calculated if no switch happened, (which shouldnt be reached anyway)

        ! in this routine i always want to include the inputted end index
        ! but only the +1 spatial orbital above semi!
        if (semi >= ende) then
            orb = nSpatOrbs + 2
            return
        end if

        ! also implement this with the new fortran 2008 routines!
        ! make the spin_change bit-rep

        orb = nSpatOrbs + 2

        do iOrb = ende, semi + 1, -1
            a = getStepvalue(ilutI, iOrb)
            b = getStepvalue(ilutJ, iOrb)
            if (a /= b) then
                orb = iOrb
                return
            end if
        end do

    end function findLastSwitch

    ! write custom add_ilut_lists
    subroutine add_guga_lists(nDets1, nDets2, list1, list2)
        integer, intent(inout) :: nDets1
        integer, intent(in) :: nDets2
        integer(n_int), intent(inout) :: list1(0:, 1:), list2(0:, 1:)

        integer :: i, min_ind, pos, abs_pos, j
        HElement_t(dp) :: tmp_mat

        ! first sort lists to use binary search
        call sort(list1(:, 1:nDets1), ilut_lt, ilut_gt)
        call sort(list2(:, 1:nDets2), ilut_lt, ilut_gt)

        abs_pos = 0
        min_ind = 1

        do i = 1, nDets2

            pos = binary_search(list1(0:nifd, min_ind:ndets1), list2(0:nifd, i))
            if (pos > 0) then
                ! try new implementation of that without the need of an extra
                ! output list

                ! need the absolute position after binary searches in only
                ! sublists
                abs_pos = abs_pos + pos

                ! when element found just update the matrix element and update
                ! the indices
                tmp_mat = extract_h_element(list1(:, abs_pos)) + &
                          extract_h_element(list2(:, i))

                call encode_matrix_element(list1(:, abs_pos), tmp_mat, 1)

                ! min_ind to search next element is then
                min_ind = min_ind + pos

            else
                ! if the entry is not in list1 i have to move down all
                ! subsequent elements after the found position and insert the
                ! new entry at the indicated absolute position
                abs_pos = abs_pos - pos

                ! is this too big to copy in one go??
                do j = nDets1, abs_pos, -1
                    list1(:, j + 1) = list1(:, j)
                end do
                ! and add the new entry
                list1(:, abs_pos) = list2(:, i)

                ! the minimum index to search from now on should not include
                ! the inserted element, since there should be no repetitions
                ! in both lists
                min_ind = min_ind - pos

                ! and i have to update the number of determinants in list1
                nDets1 = nDets1 + 1

            end if
        end do

    end subroutine add_guga_lists

    ! have to write custom x0 and x1 matrix encoding functions for the
    ! GUGA excitation integer lists, as i need an additional entry for the
    ! x1 matrix element

    function csf_purify(sd_hilbert_space, total_spin, n_el) result(csfs)
        ! function to filter out all spin-allowed states from a
        ! SD Hilbert space
        integer(n_int), intent(in) :: sd_hilbert_space(:,:)
        integer, intent(in) :: total_spin, n_el

        integer(n_int), allocatable :: csfs(:,:)
        integer :: i, cnt
        integer(n_int), allocatable :: temp_csfs(:,:)

        ! we have definitely <= sds
        allocate(temp_csfs(size(sd_hilbert_space,1),size(sd_hilbert_space,2)), &
            source = 0_n_int)

        cnt = 0

        do i = 1, size(sd_hilbert_space,2)
            if (isProperCSF_flexible(sd_hilbert_space(:,i), total_spin, n_el)) then
                cnt = cnt + 1
                temp_csfs(:,cnt) = sd_hilbert_space(:,i)
            end if
        end do

        allocate(csfs(size(sd_hilbert_space,1), cnt), source = temp_csfs(:,1:cnt))

    end function csf_purify

    pure real(dp) function get_preceeding_opposites(nJ, orb)
        integer, intent(in) :: nJ(nel), orb

        integer :: i, tgt_spin

        get_preceeding_opposites = 0.0_dp

        tgt_spin = get_spin(nJ(orb))

        do i = 1, orb - 1
            if (get_spin(nJ(i)) /= tgt_spin) then
                get_preceeding_opposites = get_preceeding_opposites + 1.0_dp
            end if
        end do

    end function get_preceeding_opposites

    subroutine write_guga_list(nunit, ilut, n_orbs)
        integer(n_int), intent(in) :: ilut(:, :)
        integer, intent(in) :: nunit
        integer, intent(in), optional :: n_orbs

        integer :: i, n_orbs_
        def_default(n_orbs_, n_orbs, nSpatorbs)

        print *, " ilut list: "
        print *, " ==========="
        do i = 1, size(ilut, 2)
            call write_det_guga(nunit, ilut(:, i), n_orbs = n_orbs_)
        end do
        print *, " ==========="

    end subroutine write_guga_list

    subroutine write_det_guga(nunit, ilut, flag, n_orbs)
        ! subroutine which prints the stepvector representation of an ilut
        integer(n_int), intent(in) :: ilut(0:)
        integer, intent(in) :: nunit
        logical, intent(in), optional :: flag
        integer, intent(in), optional :: n_orbs

        integer :: step(nSpatOrbs), i
        logical :: flag_
        integer(int_rdm) :: rdm_ind
        integer :: n_orbs_
        def_default(n_orbs_, n_orbs, nSpatorbs)
        def_default(flag_, flag, .true.)

        step = calcStepvector(ilut(0:GugaBits%len_orb))

        write(nunit, '("(")', advance='no')

        do i = 1, n_orbs_
            write(nunit, '(i3)', advance='no') step(i)
            if (i /= n_orbs_) write(nunit, '(",")', advance='no')
        end do
        write(nunit, '(")")', advance='no')

        write(nunit, '("(")', advance='no')
        do i = 1, 2
            write(nunit, "(f16.7)", advance='no') extract_matrix_element(ilut, i)
            if (i /= 2) write(nunit, "(A)", advance='no') ","
        end do

        ! if we have more entries due to RDMs, print it here
        if (ubound(ilut, dim = 1) == GugaBits%len_tot) then
            if (tRDMonfly) then
                rdm_ind = extract_rdm_ind(ilut)
                write(nunit, "(A,i8)", advance='no') ") ", getDeltaB(ilut)
                write(nunit, "(A,3i8,A)", advance='no') " | ( ", &
                    extract_excit_lvl_rdm(rdm_ind), &
                    extract_excit_type_rdm(rdm_ind), &
                    int(iand(rdm_ind, rdm_ind_bitmask)), ' ) '
                write(nunit, "(f16.7)", advance='no') &
                    extract_stochastic_rdm_x0(GugaBits, ilut)
                if (flag_) then
                    write(nunit, "(f16.7)", advance='yes') &
                        extract_stochastic_rdm_x1(GugaBits, ilut)
                else
                    write(nunit, "(f16.7)", advance='no') &
                        extract_stochastic_rdm_x1(GugaBits, ilut)
                end if
            else
                if (flag_) then
                    write(nunit, "(A,i8)", advance='yes') ") ", getDeltaB(ilut)
                else
                    write(nunit, "(A,i8)", advance='no') ") ", getDeltaB(ilut)
                end if
            end if
        else
            if (flag_) then
                write(nunit, "(A,i8)", advance = 'yes') ") ", getDeltaB(ilut)
            else
                write(nunit, "(A)", advance='no') ") "
            end if
        end if

    end subroutine write_det_guga

    subroutine encode_matrix_element_real(ilut, mat_ele, mat_type)
        ! encodes the x0 or x1 matrix element needed during the excitation
        ! creation.
        ! mat_ele   ... x0 or x1 matrix element
        ! mat_type  ... 1...x0, 2...x1
        integer(n_int), intent(inout) :: ilut(0:GugaBits%len_tot)
        real(dp), intent(in) :: mat_ele
        integer, intent(in) :: mat_type
        character(*), parameter :: this_routine = "encode_matrix_element_real"

        integer(n_int) :: mat_int ! integer version of real

        ASSERT(mat_type == 1 .or. mat_type == 2)

        mat_int = transfer(mat_ele, mat_int)

        ilut(GugaBits%len_orb + mat_type) = mat_int

    end subroutine encode_matrix_element_real

#ifdef CMPLX_
    subroutine encode_matrix_element_cmplx(ilut, mat_ele, mat_type)
        ! this is specific for complex matrix elements.. here
        ! i can use the two storage slots for x0 and x1 to encode
        ! both the real and imaginary parts of the Hamiltonian matrix elements
        integer(n_int), intent(inout) :: ilut(0:GugaBits%len_tot)
        complex(dp), intent(in) :: mat_ele
        integer, intent(in) :: mat_type
        character(*), parameter :: this_routine = "encode_matrix_element_cmplx"

        integer(n_int) :: mat_int

        ! here I have to assert that mat_ind is 1! otherwise smth went wrong
        ASSERT(mat_type == 1)
        ! also make sure the x1 mat ele was nullified.
!         ASSERT(near_zero(extract_matrix_element(ilut, 2)))

        mat_int = transfer(real(real(mat_ele), dp), mat_int)
        ilut(GugaBits%ind_x0) = mat_int

        mat_int = transfer(real(aimag(mat_ele), dp), mat_int)
        ilut(GugaBits%ind_x1) = mat_int

    end subroutine encode_matrix_element_cmplx
#endif

    function extract_matrix_element(ilut, mat_type) result(mat_ele)
        ! function to extract matrix element of a GUGA ilut
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        integer, intent(in) :: mat_type
        real(dp) :: mat_ele
        character(*), parameter :: this_routine = "extract_matrix_element"

        ASSERT(mat_type == 1 .or. mat_type == 2)

        mat_ele = transfer(ilut(GugaBits%len_orb + mat_type), mat_ele)

    end function extract_matrix_element

    subroutine update_matrix_element_real(ilut, mat_ele, mat_type)
        ! function to update already encoded matrix element multiplicative
        integer(n_int), intent(inout) :: ilut(0:GugaBits%len_tot)
        real(dp), intent(in) :: mat_ele
        integer, intent(in) :: mat_type
        character(*), parameter :: this_routine = "update_matrix_element_real"

        integer(n_int) :: mat_int
        real(dp) :: temp_ele

        ASSERT(mat_type == 1 .or. mat_type == 2)

        temp_ele = transfer(ilut(GugaBits%len_orb + mat_type), temp_ele)

        mat_int = transfer(temp_ele * mat_ele, mat_int)

        ilut(GugaBits%len_orb + mat_type) = mat_int

    end subroutine update_matrix_element_real

#ifdef CMPLX_
    subroutine update_matrix_element_cmplx(ilut, mat_ele, mat_type)
        ! specific function if we need to update with a complex integral
        integer(n_int), intent(inout) :: ilut(0:GugaBits%len_tot)
        complex(dp), intent(in) :: mat_ele
        integer, intent(in) :: mat_type
        character(*), parameter :: this_routine = "update_matrix_element_cmplx"

        integer(n_int) :: mat_int
        real(dp) :: temp_ele

        ! here i have to again assert that mat_type == 1, otherwise
        ! smth went wrong..
        ASSERT(mat_type == 1)

        ! just for checking see if the x2 is nullified always..
        ASSERT(near_zero(extract_matrix_element(ilut, 2)))

        ! and now we want to store the real and imag part in x0 and x1..
        temp_ele = transfer(ilut(GugaBits%ind_x0), temp_ele)
        mat_int = transfer(temp_ele * real(real(mat_ele), dp), mat_int)

        ilut(GugaBits%ind_x0) = mat_int

        mat_int = transfer(temp_ele * real(aimag(mat_ele), dp), mat_int)
        ilut(GugaBits%ind_x1) = mat_int

    end subroutine update_matrix_element_cmplx
#endif

    function count_beta_orbs_ij(csf_i, i, j) result(nOpen)
        ! function to count the number of 1s in a CSF det between spatial
        ! orbitals i and j
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: i, j
        integer :: nOpen
        character(*), parameter :: this_routine = "count_beta_orbs_ij"

        integer :: k

        ASSERT(i > 0 .and. i <= nSpatOrbs)
        ASSERT(j > 0 .and. j <= nSpatOrbs)

        nOpen = 0

        ! quick and dirty fix to deal with the excitation range mask probs:
        ! do i always call that for the current det the excitation is
        ! calculated for?  i think so..
        do k = i, j
            if (csf_i%stepvector(k) == 1) then
                nOpen = nOpen + 1
            end if
        end do
    end function count_beta_orbs_ij

    function count_alpha_orbs_ij(csf_i, i, j) result(nOpen)
        ! function to count the number of 2s in a CSF det between spatial
        ! orbitals i and j
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: i, j
        integer :: nOpen
        character(*), parameter :: this_routine = "count_alpha_orbs_ij"

        integer :: k

        ASSERT(i > 0 .and. i <= nSpatOrbs)
        ASSERT(j > 0 .and. j <= nSpatOrbs)

        nOpen = 0

        ! quick fix for now to see if thats the problem: loop and check!
        do k = i, j
            if (csf_i%stepvector(k) == 2) then
                nOpen = nOpen + 1
            end if
        end do
    end function count_alpha_orbs_ij

    function count_open_orbs_ij(csf_i, i, j, L) result(nOpen)
        ! function to calculate the number of open orbitals between spatial
        ! orbitals i and j in ilut. i and j have to be given ordered i<j
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: i, j
        integer(n_int), intent(in), optional :: L(0:GugaBits%len_orb)
        integer :: nOpen
        character(*), parameter :: this_routine = "count_open_orbs_ij"

        logical :: flag
        integer :: k

        ASSERT(i > 0 .and. i <= nSpatOrbs)
        ASSERT(j > 0 .and. j <= nSpatOrbs)
        ! scrap this assert and change in that way to output 0 if the indices
        ! dont fit or are reversed. to deal with to short overlap ranges

        nOpen = 0

        ! also here a quick fix do deal with excitrangemask probs:

        ! if the ilut input is present use it otherwise just look at the
        ! stepvector
        if (present(L)) then
            do k = i, j
                flag = isOne(L, k)
                if (flag .or. isTwo(L, k)) then
                    nOpen = nOpen + 1
                end if
            end do
        else
            do k = i, j
                if (csf_i%stepvector(k) == 1 .or. csf_i%stepvector(k) == 2) then
                    nOpen = nOpen + 1
                end if
            end do
        end if

    end function count_open_orbs_ij

    function getExcitationRangeMask(i, j) result(mask)
        ! function to create an integer corresponding to a mask were every
        ! bits between the spin orbitals corresponding to spatial orbitals
        ! i and j are set to 1 and 0 else. used to access the excitation
        ! range for GUGA excitation calculations
        integer, intent(in) :: i, j
        integer(n_int) :: mask
        character(*), parameter :: this_routine = "getExcitationRangeMask"

        integer(n_int) :: k
        integer(n_int) :: tmp_i, tmp_j

        ASSERT(i > 0 .and. i <= nSpatOrbs)
        ASSERT(j > 0 .and. j <= nSpatOrbs)
        ASSERT(j > i)

        tmp_i = int(i, n_int)
        tmp_j = int(j, n_int)

        ! not quite sure about LMSB or RMSB... todo
        mask = 0_n_int

        do k = 2_n_int * tmp_i - 1_n_int, 2_n_int * tmp_j ! convert to spin orbitals
            mask = mask + 2_n_int**(k - 1_n_int)
        end do

    end function getExcitationRangeMask

    pure subroutine setDeltaB(deltaB, ilut)
        ! routine to encode the deltaB value in a given CSF in ilut bit
        ! representation, by using the newly defined flags:
        ! flag_deltaB_sign   ... 7
        ! flag_deltaB_single ... 5
        ! and if necessary
        ! flag_deltaB_double ... 6
        integer, intent(in) :: deltaB
        integer(n_int), intent(inout) :: ilut(0:GugaBits%len_tot)

        ! should no just be:
        ilut(GugaBits%ind_b) = deltaB

    end subroutine setDeltaB

    function getDeltaB(ilut) result(deltaB)
        ! function to get the deltaB value encoded in the flag-byte in ilut
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        integer :: deltaB

        ! check if flags are correctly set

        ! and this should now jsut be:
        deltaB = int(ilut(GugaBits%ind_b))

    end function getDeltaB

    function extract_h_element(ilutG) result(HElement)
        integer(n_int), intent(in) :: ilutG(0:GugaBits%len_tot)
        HElement_t(dp) :: HElement

#ifdef CMPLX_
        HElement = cmplx(extract_matrix_element(ilutG, 1), &
                         extract_matrix_element(ilutG, 2), kind=dp)
#else
        HElement = extract_matrix_element(ilutG, 1)
#endif

    end function extract_h_element

    subroutine convert_ilut_toNECI(ilutG, ilutN, HElement)
        integer(n_int), intent(in) :: ilutG(0:GugaBits%len_tot)
        integer(n_int), intent(out) :: ilutN(0:niftot)
        HElement_t(dp), intent(out), optional :: HElement
        character(*), parameter :: this_routine = "convert_ilut_toNECI"

        ASSERT(isProperCSF_ilut(ilutG))

        ilutN = 0_n_int
        ! i think i just need to copy over the det part again
        ilutN(0:GugaBits%len_orb) = ilutG(0:GugaBits%len_orb)

        ! and then extract the matrix element
        ! here i need to check what type of matrix element is necessary
        ! dependent on which type of compilation,
        ! extract_matrix_element always gives a real(dp)!
        if (present(HElement)) then
            HElement = extract_h_element(ilutG)
        end if

        if (tFillingStochRDMonfly) then
            ! in this case I need to transfer the rdm_ind and x0,x1 info to
            ! the 'neci ilut'
            call transfer_stochastic_rdm_info(ilutG, ilutN)
        end if

    end subroutine convert_ilut_toNECI

    pure subroutine transfer_stochastic_rdm_info(ilutG, ilutN, &
                                                 BitIndex_from, BitIndex_to)
        ! I need to keep this general with BitIndex unfortunately because
        ! i have to perform this on ilut with niftot and on parent arrays..
        integer(n_int), intent(in) :: ilutG(0:)
        integer(n_int), intent(inout) :: ilutN(0:)
        type(BitRep_t), intent(in), optional :: BitIndex_from, BitIndex_to

        type(BitRep_t) :: from, to

        def_default(from, BitIndex_from, GugaBits)
        def_default(to, BitIndex_to, IlutBits)

        ! here i now I get a GUGA ilut, but could be that I have to
        ! transfer to Parent array? I dont think so..

        call encode_stochastic_rdm_info(to, ilutN, &
                                        rdm_ind=extract_stochastic_rdm_ind(from, ilutG), &
                                        x0=extract_stochastic_rdm_x0(from, ilutG), &
                                        x1=extract_stochastic_rdm_x1(from, ilutG))

    end subroutine transfer_stochastic_rdm_info

    subroutine convert_ilut_toGUGA(ilutN, ilutG, HElement, delta_b)
        integer(n_int), intent(in) :: ilutN(0:niftot)
        integer(n_int), intent(out) :: ilutG(0:GugaBits%len_tot)
        HElement_t(dp), intent(in), optional :: HElement
        integer, intent(in), optional :: delta_b

        ilutG = 0_n_int

        ! need only the det part essentially..
        ilutG(0:GugaBits%len_orb) = ilutN(0:GugaBits%len_orb)

        if (present(HElement)) then
            call encode_matrix_element(ilutG, 0.0_dp, 2)
            call encode_matrix_element(ilutG, HElement, 1)
        else
            ! and set matrix elements to 1 and delta b to 0
            call encode_matrix_element(ilutG, 1.0_dp, 1)
            call encode_matrix_element(ilutG, 0.0_dp, 2)
        end if

        if (present(delta_b)) then
            call setDeltaB(delta_b, ilutG)
        else
            call setDeltaB(0, ilutG)
        end if

    end subroutine convert_ilut_toGUGA

    function isProperCSF_nI(nI) result(flag)
        ! function to check if provided CSF in nI(nEl) format is a proper CSF
        integer, intent(in) :: nI(nEl)
        logical :: flag

        real(dp) :: bVector(nEl)
        integer :: calcS

        flag = .true.

        bVector = calcB_vector_nI(nI)

        ! check if the b value drops below zero
        if (any(bVector < 0.0_dp)) flag = .false.

        ! check if total spin is same as input
        calcS = abs(sum(get_spin_pn(nI(:))))

        if (calcS /= STOT) flag = .false.

    end function isProperCSF_nI

    pure function isProperCSF_b(ilut) result(flag)
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        logical :: flag

        flag = all(calcB_vector_int(ilut(0:GugaBits%len_orb)) >= 0)
    end function isProperCSF_b

    pure function isProperCSF_flexible(ilut, spin, num_el) result(flag)
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        integer, intent(in) :: spin, num_el
        logical :: flag

        flag = (all(calcB_vector_int(ilut(0:GugaBits%len_orb)) > 0) &
            .and. (abs(return_ms(ilut, num_el)) == spin) &
            .and. (int(sum(calcOcc_vector_ilut(ilut(0:GugaBits%len_orb)))) == num_el))

    end function isProperCSF_flexible

    function isProperCSF_sys(ilut, sysFlag, t_print_in) result(flag)
        ! function to check if provided CSF in ilut format is a proper CSF
        ! checks b vector positivity and is total S is correct
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        logical, intent(in):: sysFlag
        logical, intent(in), optional :: t_print_in
        logical :: flag

        logical :: t_print

        if (present(t_print_in)) then
            t_print = t_print_in
        else
            t_print = .false.
        end if

        flag = .true.

        ! check if b value drops below zero
        if (any(calcB_vector_int(ilut(0:GugaBits%len_orb)) < 0)) flag = .false.

        ! if system flag is also given as input also check if the CSF fits
        ! concerning total S and the number of electrons
        if (sysFlag) then
            if (abs(return_ms(ilut)) /= STOT) then
                if (t_print) then
                    print *, "CSF does not have correct total spin!:"
                    call write_det_guga(6, ilut)
                    print *, "System S: ", STOT
                    print *, "CSF S: ", abs(return_ms(ilut))
                end if
                flag = .false.
            end if

            if (int(sum(calcOcc_vector_ilut(ilut(0:GugaBits%len_orb)))) /= nEl) then
                if (t_print) then
                    print *, "CSF does not have right number of electrons!:"
                    call write_det_guga(6, ilut)
                    print *, "System electrons: ", nEl
                    print *, "CSF electrons: ", &
                        int(sum(calcOcc_vector_ilut(ilut(0:GugaBits%len_orb))))
                end if
                flag = .false.
            end if
        end if

    end function isProperCSF_sys

    pure function calcB_vector_nI(nI) result(bVector)
        ! function to calculate the bVector from a CSF given in nI
        ! representation. Gives b vector also of length nEl.
        ! not yet quite sure if i should output b as integers or real
        integer, intent(in) :: nI(nEl)
        real(dp) :: bVector(nEl), bValue

        integer :: iOrb, inc

        ! init
        iOrb = 1
        bVector = 0.0_dp
        bValue = 0.0_dp
        inc = 0

        do while (iOrb <= nEl)
            if (isDouble(nI, iOrb)) then
                bVector(iOrb) = bValue
                bVector(iOrb + 1) = bValue
                inc = 2

            else
                inc = 1
                if is_alpha(nI(iOrb)) then
                bValue = bValue - 1.0_dp

            else
                bValue = bValue + 1.0_dp

            end if
            bVector(iOrb) = bValue

            end if

            iOrb = iOrb + inc

        end do

    end function calcB_vector_nI

    pure function isDouble(nI, iOrb) result(flag)
        ! returns a logical .true. if spinorbital iOrb is doubly occupied
        ! and .false. elsewise.
        integer, intent(in) :: nI(nEl), iOrb
        logical :: flag

        integer :: pair

        flag = .false.
        pair = ab_pair(nI(iOrb))

        ! ask simon about the ordering in nI. If its not always ordered I
        ! have to to a quicksearch for pair in nI. If its ordered I can
        ! just check adjacent entries in nI

        ! assume ordered for now!
        if (iOrb == 1) then
            ! special conditions if iOrb is 1
            flag = (nI(2) == pair)
        else if (iOrb == nEl) then
            flag = (nI(nEl - 1) == pair)
        else
            flag = ((nI(iOrb - 1) == pair) .or. (nI(iOrb + 1) == pair))
        end if

    end function isDouble

    pure function calcStepvector(ilut) result(stepVector)
        ! function to calculate stepvector of length nReps, corresponding
        ! to the ilut bit-representation, if each stepvalue is needed
        ! often within a function.
        ! there is probably a very efficient way of programming that!
        !TODO ask simon for improvements.
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_orb)
        integer :: stepVector(nSpatOrbs)
        integer :: iOrb

        do iOrb = 1, nSpatOrbs
            stepVector(iOrb) = getStepvalue(ilut, iOrb)
        end do

    end function calcStepvector

    pure function calcOcc_vector_ilut(ilut) result(occVector)
        ! probably more efficiently implemented by simon already...
        ! but for now do it in this stupid way todo -> ask simon
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_orb)
        real(dp) :: occVector(nSpatOrbs)

        integer :: iOrb

        do iOrb = 1, nSpatOrbs
            occVector(iOrb) = getSpatialOccupation(ilut, iOrb)
        end do

    end function calcOcc_vector_ilut

    pure function calcOcc_vector_int(ilut) result(occVector)
        ! function which gives the occupation vector in integer form
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_orb)
        integer :: occVector(nSpatOrbs)

        integer :: i

        do i = 1, nSpatOrbs
            occVector(i) = int(getSpatialOccupation(ilut, i))
        end do

    end function calcOcc_vector_int

    pure function calcB_vector_ilut(ilut) result(bVector)
        ! function to calculate bVector of length (nBasis) for a given
        ! CSF bitRep ilut of length (2*nBasis), save this b-vector
        ! in the nI array, normally used to store occupied orbitals
        ! update: changed convention to not use nI for bVector but also for
        ! occupied orbitals as in normal NECI implementation
        ! but nethertheless need a b-vector of lenght of spatial orbitals
        ! for excitaiton and matrix element calculation
        ! UPDATE: change b vector calculation to update b vector directly on
        ! the spot of the corresponding stepvector, eg:
        ! d = 0 1 2 3
        ! b = 0 1 0 0
        ! because otherwise i actually only needed the bvector of one orbital
        ! ahead. and never the first entry otherwise. and that would cause
        ! problems when accessing b vector at the end of an excitaiton if
        ! that is the last orbital..
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_orb)
        real(dp) :: bVector(nSpatOrbs)        ! b-vector stored in bVector
        integer :: i
        real(dp) :: bValue

        ! loop over CSF entries and increment, decrement bValue
        ! accordingly, have to correctly access ilut entries and
        ! check if they are 0,1,2 or 3... -> how to for multiple
        ! integers?
        bVector = 0.0_dp
        bValue = 0.0_dp

        do i = 1, nSpatOrbs

            if (isOne(ilut, i)) then
                bValue = bValue + 1.0_dp
            else if (isTwo(ilut, i)) then
                bValue = bValue - 1.0_dp
            end if
            ! define bvalue to always only get updated for the next
            ! UPDATE: changed definition to update on the spot.
            bVector(i) = bValue
        end do

    end function calcB_vector_ilut

    pure function calcB_vector_int(ilut) result(bVector)
        ! function to calculate the bvector in integer form
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_orb)
        integer :: bVector(nSpatOrbs)

        integer :: i, bValue

        bVector = 0
        bValue = 0

        do i = 1, nSpatOrbs

            if (isOne(ilut, i)) then
                bValue = bValue + 1

            else if (isTwo(ilut, i)) then
                bValue = bValue - 1

            end if

            bVector(i) = bValue

        end do

    end function calcB_vector_int
    pure subroutine EncodeBitDet_guga(nI, ilut)
        ! special function to encode bit dets for the use in the guga
        ! excitation generation
        integer, intent(in) :: nI(nEl)
        integer(n_int), intent(out) :: ilut(0:GugaBits%len_tot)

        integer :: i, pos

        ilut = 0_n_int

        do i = 1, nEl
            pos = (nI(i) - 1) / bits_n_int
            ilut(pos) = ibset(ilut(pos), mod(nI(i) - 1, bits_n_int))
        end do

    end subroutine EncodeBitDet_guga

    pure function getSpatialOccupation(iLut, s) result(nOcc)

        integer(n_int), intent(in) :: ilut(0:GugaBits%len_orb)
        integer, intent(in) :: s
        real(dp) :: nOcc

        if (isZero(ilut, s)) then
            nOcc = 0.0_dp

        else if (isThree(ilut, s)) then
            nOcc = 2.0_dp

        else
            nOcc = 1.0_dp
        end if

    end function getSpatialOccupation

    pure function convert_guga_to_ni(csf, siz) result(nI)
        ! function to make it easier for me to input a csf in my used notation
        ! to a nI NECI array..
        integer, intent(in) :: siz
        integer, intent(in) :: csf(siz)
        integer :: nI(nel)

        integer :: i, cnt_orbs, cnt_ind

        cnt_orbs = 0
        cnt_ind = 0

        do i = 1, siz

            select case (csf(i))

            case (0)
                ! nothing to do actually except update the

            case (1)
                ! beta orbital
                cnt_ind = cnt_ind + 1

                nI(cnt_ind) = cnt_orbs + 1

            case (2)
                ! alpha orbs

                cnt_ind = cnt_ind + 1

                nI(cnt_ind) = cnt_orbs + 2

            case (3)
                ! doubly occupied
                cnt_ind = cnt_ind + 1
                nI(cnt_ind) = cnt_orbs + 1
                cnt_ind = cnt_ind + 1
                nI(cnt_ind) = cnt_orbs + 2

            end select

            ! update orbitals for every csf entry
            cnt_orbs = cnt_orbs + 2

        end do

    end function convert_guga_to_ni

    pure subroutine calc_csf_i(ilut, step_vector, b_vector, occ_vector)
        ! routine to calculate the csf information for specific outputted
        ! information
        integer(n_int), intent(in) :: ilut(0:niftot)
        integer, intent(out) :: step_vector(nSpatOrbs), b_vector(nSpatOrbs)
        real(dp), intent(out) :: occ_vector(nSpatOrbs)

        integer :: b_int, i, step
        ! copy the stuff from below.. when do i want to allocate the objects?
        ! hm..
        step_vector = 0
        b_vector = 0
        occ_vector = 0.0_dp

        b_int = 0

        do i = 1, nSpatOrbs
            step = getStepvalue(ilut, i)
            step_vector(i) = step
            select case (step)
            case (1)
                occ_vector(i) = 1.0_dp
                b_int = b_int + 1
            case (2)
                occ_vector(i) = 1.0_dp
                b_int = b_int - 1
            case (3)
                occ_vector(i) = 2.0_dp
            end select

            b_vector(i) = b_int
        end do
    end subroutine calc_csf_i

    pure function construct_CSF_Info_t(ilut) result(csf_i)
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        type(CSF_Info_t) :: csf_i
        call new_CSF_Info_t(nSpatOrbs, csf_i)
        call fill_csf_i(ilut, csf_i)
    end function

    pure subroutine new_CSF_Info_t(n_spat_orbs, csf_i)
        integer, intent(in) :: n_spat_orbs
        type(CSF_Info_t), intent(out) :: csf_i
        allocate(csf_i%stepvector(n_spat_orbs), &
                 csf_i%B_real(n_spat_orbs), &
                 csf_i%Occ_real(n_spat_orbs), &
                 csf_i%B_int(n_spat_orbs), &
                 csf_i%Occ_int(n_spat_orbs), &
                 csf_i%cum_list(n_spat_orbs))
    end subroutine

    pure subroutine fill_csf_i(ilut, csf_i)
        ! routine which sets up all the additional csf information, like
        ! stepvector, b vector, occupation etc. in various formats in one
        ! place
        ! and combine all the necessary calcs. into one loop instead of
        ! the seperate ones..
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        type(CSF_Info_t), intent(inout) :: csf_i
        debug_function_name("fill_csf_i")

        integer :: i, step, b_int
        real(dp) :: b_real, cum_sum

        ASSERT(isProperCSF_ilut(ilut))
        ASSERT(allocated(csf_i%stepvector))
        ASSERT(allocated(csf_i%B_real))
        ASSERT(allocated(csf_i%Occ_real))
        ASSERT(allocated(csf_i%B_int))
        ASSERT(allocated(csf_i%Occ_int))

        csf_i%stepvector = 0
        csf_i%B_real = 0.0_dp
        csf_i%Occ_real = 0.0_dp
        csf_i%B_int = 0
        csf_i%Occ_int = 0

        b_real = 0.0_dp
        b_int = 0


        ! TODO(@Oskar): Use these functions instead
        ! currentB_ilut = calcB_vector_ilut(ilut)
        ! currentOcc_ilut = calcOcc_vector_ilut(ilut)
        ! currentOcc_int = calcOcc_vector_int(ilut)
        ! current_stepvector = calcStepVector(ilut)
        ! currentB_int = calcB_vector_int(ilut)

        ! also create a fake cum-list of the non-doubly occupied orbitals
        csf_i%cum_list = 0.0_dp
        cum_sum = 0.0_dp

        do i = 1, nSpatOrbs

            step = getStepvalue(ilut, i)

            csf_i%stepvector(i) = step

            select case (step)

            case (0)

                csf_i%Occ_real(i) = 0.0_dp
                csf_i%Occ_int(i) = 0

                cum_sum = cum_sum + 1.0_dp

            case (1)

                csf_i%Occ_real(i) = 1.0_dp
                csf_i%Occ_int(i) = 1

                b_real = b_real + 1.0_dp
                b_int = b_int + 1

                cum_sum = cum_sum + 1.0_dp

            case (2)

                csf_i%Occ_real(i) = 1.0_dp
                csf_i%Occ_int(i) = 1

                b_real = b_real - 1.0_dp
                b_int = b_int - 1

                cum_sum = cum_sum + 1.0_dp

            case (3)

                csf_i%Occ_real(i) = 2.0_dp
                csf_i%Occ_int(i) = 2

            end select

            csf_i%B_real(i) = b_real
            csf_i%B_int(i) = b_int

            csf_i%cum_list(i) = cum_sum
        end do

    end subroutine fill_csf_i

    pure function is_compatible(ilut, csf_i) result(res)
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        type(CSF_Info_t), intent(in) :: csf_i
        logical :: res
        ! TODO(@Oskar): Implement
        unused_var(ilut); unused_var(csf_i)
        res = .true.
    end function

    pure subroutine encode_stochastic_rdm_info(BitIndex, ilut, rdm_ind, x0, x1)
        ! make these function general by also providing the
        ! bit-rep index data-structure!
        type(BitRep_t), intent(in) :: BitIndex
        integer(n_int), intent(inout) :: ilut(0:BitIndex%len_tot)
        integer(int_rdm), intent(in) :: rdm_ind
        real(dp), intent(in) :: x0, x1

        ! i need to be sure that int_rdm and n_int are of the same
        ! size otherwise this breaks..
        call encode_stochastic_rdm_ind(BitIndex, ilut, rdm_ind)
        call encode_stochastic_rdm_x0(BitIndex, ilut, x0)
        call encode_stochastic_rdm_x1(BitIndex, ilut, x1)

    end subroutine encode_stochastic_rdm_info

    pure subroutine extract_stochastic_rdm_info(BitIndex, ilut, rdm_ind, x0, x1)
        type(BitRep_t), intent(in) :: BitIndex
        integer(n_int), intent(in) :: ilut(0:BitIndex%len_tot)
        integer(int_rdm), intent(out) :: rdm_ind
        real(dp), intent(out) :: x0, x1

        rdm_ind = extract_stochastic_rdm_ind(BitIndex, ilut)
        x0 = extract_stochastic_rdm_x0(BitIndex, ilut)
        x1 = extract_stochastic_rdm_x1(BitIndex, ilut)

    end subroutine extract_stochastic_rdm_info

    pure subroutine encode_stochastic_rdm_ind(BitIndex, ilut, rdm_ind)
        type(BitRep_t), intent(in) :: BitIndex
        integer(n_int), intent(inout) :: ilut(0:BitIndex%len_tot)
        integer(int_rdm), intent(in) :: rdm_ind

        ilut(BitIndex%ind_rdm_ind) = rdm_ind

    end subroutine encode_stochastic_rdm_ind

    pure function extract_stochastic_rdm_ind(BitIndex, ilut) result(rdm_ind)
        type(BitRep_t), intent(in) :: BitIndex
        integer(n_int), intent(in) :: ilut(0:BitIndex%len_tot)
        integer(int_rdm) :: rdm_ind

        rdm_ind = ilut(BitIndex%ind_rdm_ind)

    end function extract_stochastic_rdm_ind

    pure subroutine encode_stochastic_rdm_x0(BitIndex, ilut, x0)
        type(BitRep_t), intent(in) :: BitIndex
        integer(n_int), intent(inout) :: ilut(0:BitIndex%len_tot)
        real(dp), intent(in) :: x0

        integer(n_int) :: x0_int

        x0_int = transfer(x0, x0_int)

        ilut(BitIndex%ind_rdm_x0) = x0_int

    end subroutine encode_stochastic_rdm_x0

    pure function extract_stochastic_rdm_x0(BitIndex, ilut) result(x0)
        type(BitRep_t), intent(in) :: BitIndex
        integer(n_int), intent(in) :: ilut(0:BitIndex%len_tot)
        real(dp) :: x0

        x0 = transfer(ilut(BitIndex%ind_rdm_x0), x0)

    end function extract_stochastic_rdm_x0

    pure subroutine encode_stochastic_rdm_x1(BitIndex, ilut, x1)
        type(BitRep_t), intent(in) :: BitIndex
        integer(n_int), intent(inout) :: ilut(0:BitIndex%len_tot)
        real(dp), intent(in) :: x1

        integer(n_int) :: x1_int

        x1_int = transfer(x1, x1_int)

        ilut(BitIndex%ind_rdm_x1) = x1_int

    end subroutine encode_stochastic_rdm_x1

    pure function extract_stochastic_rdm_x1(BitIndex, ilut) result(x1)
        type(BitRep_t), intent(in) :: BitIndex
        integer(n_int), intent(in) :: ilut(0:BitIndex%len_tot)
        real(dp) :: x1

        x1 = transfer(ilut(BitIndex%ind_rdm_x1), x1)

    end function extract_stochastic_rdm_x1

    pure subroutine extract_1_rdm_ind(rdm_ind, i, a, excit_lvl, excit_typ)
        ! the converstion routine between the combined and explicit rdm
        ! indices for the 1-RDM
        integer(int_rdm), intent(in) :: rdm_ind
        integer, intent(out) :: i, a
        integer, intent(out), optional :: excit_lvl, excit_typ

        integer(int_rdm) :: rdm_ind_

        ! if we also want to use the top 7 bits of rdm_ind for information
        ! of the excit-lvl and type we have to 0 them out before
        ! extracting the indices

        rdm_ind_ = iand(rdm_ind, rdm_ind_bitmask)

        a = int(mod(rdm_ind_ - 1, nSpatOrbs) + 1)
        i = int((rdm_ind_ - 1) / nSpatOrbs + 1)

        if (present(excit_lvl)) then
            excit_lvl = extract_excit_lvl_rdm(rdm_ind)
        end if

        if (present(excit_typ)) then
            excit_typ = extract_excit_type_rdm(rdm_ind)
        end if

    end subroutine extract_1_rdm_ind

    pure function contract_1_rdm_ind(i, a, excit_lvl, excit_typ) result(rdm_ind)
        ! the inverse function of the routine above, to give the combined
        ! rdm index of two explicit ones
        integer, intent(in) :: i, a
        integer, intent(in), optional :: excit_lvl, excit_typ
        integer(int_rdm) :: rdm_ind

        rdm_ind = nSpatOrbs * (i - 1) + a

        if (present(excit_lvl)) then
            call encode_excit_lvl_rdm(rdm_ind, excit_lvl)
        end if

        if (present(excit_typ)) then
            call encode_excit_typ_rdm(rdm_ind, excit_typ)
        end if

    end function contract_1_rdm_ind

    pure function extract_excit_type_rdm(rdm_ind) result(excit_typ)
        integer(int_rdm), intent(in) :: rdm_ind
        integer :: excit_typ

        excit_typ = int(ibits(rdm_ind, pos_excit_type_bits, n_excit_type_bits))

    end function extract_excit_type_rdm

    pure function extract_excit_lvl_rdm(rdm_ind) result(excit_lvl)
        integer(int_rdm), intent(in) :: rdm_ind
        integer :: excit_lvl

        excit_lvl = int(ibits(rdm_ind, pos_excit_lvl_bits, n_excit_lvl_bits))

    end function extract_excit_lvl_rdm

    pure subroutine encode_excit_typ_rdm(rdm_ind, excit_typ)
        integer(int_rdm), intent(inout) :: rdm_ind
        integer, intent(in) :: excit_typ

        call mvbits(int(excit_typ, int_rdm), 0, n_excit_type_bits, &
                    rdm_ind, pos_excit_type_bits)

    end subroutine encode_excit_typ_rdm


    pure subroutine encode_excit_lvl_rdm(rdm_ind, excit_lvl)
        integer(int_rdm), intent(inout) :: rdm_ind
        integer, intent(in) :: excit_lvl

        ! i need to mv the bit-rep of excit_lvl to the corresponding
        ! position in rdm_ind
        call mvbits(int(excit_lvl, int_rdm), 0, n_excit_lvl_bits, &
                    rdm_ind, pos_excit_lvl_bits)

    end subroutine encode_excit_lvl_rdm

    pure function contract_2_rdm_ind(i, j, k, l, excit_lvl, excit_typ) result(ijkl)
        ! since I only ever have spatial orbitals in the GUGA-RDM make
        ! the definition of the RDM-index combination differently
        integer, intent(in) :: i, j, k, l
        integer, intent(in), optional :: excit_lvl, excit_typ
        integer(int_rdm) :: ijkl

        integer(int_rdm) :: ij, kl

        ij = contract_1_rdm_ind(i, j)
        kl = contract_1_rdm_ind(k, l)

        ijkl = (ij - 1) * (nSpatOrbs**2) + kl

        if (present(excit_lvl)) then
            call encode_excit_lvl_rdm(ijkl, excit_lvl)
        end if

        if (present(excit_typ)) then
            call encode_excit_typ_rdm(ijkl, excit_typ)
        end if

    end function contract_2_rdm_ind

    pure subroutine extract_2_rdm_ind(ijkl, i, j, k, l, &
                                      ij_out, kl_out, excit_lvl, excit_typ)
        ! the inverse routine of the function above.
        ! it is actually practical to have ij and kl also available at
        ! times, since it can identify diagonal entries of the two-RDM
        integer(int_rdm), intent(in) :: ijkl
        integer, intent(out) :: i, j, k, l
        integer(int_rdm), intent(out), optional :: ij_out, kl_out
        integer, intent(out), optional :: excit_lvl, excit_typ

        integer(int_rdm) :: ij, kl

        integer(int_rdm) :: ijkl_

        ijkl_ = iand(ijkl, rdm_ind_bitmask)

        kl = mod(ijkl_ - 1, int(nSpatOrbs, int_rdm)**2) + 1
        ij = (ijkl_ - kl) / (nSpatOrbs**2) + 1

        call extract_1_rdm_ind(ij, i, j)
        call extract_1_rdm_ind(kl, k, l)

        if (present(ij_out)) ij_out = ij
        if (present(kl_out)) kl_out = kl

        if (present(excit_lvl)) then
            excit_lvl = extract_excit_lvl_rdm(ijkl)
        end if
        if (present(excit_typ)) then
            excit_typ = extract_excit_type_rdm(ijkl)
        end if

    end subroutine extract_2_rdm_ind

    pure function extract_rdm_ind(ilutG) result(rdm_ind)
        integer(n_int), intent(in) :: ilutG(0:GugaBits%len_tot)
        integer(int_rdm) :: rdm_ind

        rdm_ind = ilutG(GugaBits%ind_rdm_ind)

    end function extract_rdm_ind

    pure subroutine encode_rdm_ind(ilutG, rdm_ind)
        integer(n_int), intent(inout) :: ilutG(0:GugaBits%len_tot)
        integer(int_rdm), intent(in) :: rdm_ind

        ilutG(GugaBits%ind_rdm_ind) = rdm_ind

    end subroutine encode_rdm_ind
    function encode_excit_info_vec(typ, inds) result(excit_info_int)
        ! function to encode the minimal information of an excit-info
        ! object into a single 64bit integer. used in the PCHB excitation
        ! generation.
        debug_function_name("encode_excit_info_vec")
        integer, intent(in) :: typ, inds(4)
        integer(int64) :: excit_info_int

#ifdef DEBUG_
        select case(typ)
        case( excit_type%single_overlap_L_to_R)
        case( excit_type%single_overlap_R_to_L )
        case( excit_type%double_lowering )
        case( excit_type%double_raising )
        case( excit_type%double_L_to_R_to_L)
        case( excit_type%double_R_to_L_to_R )
        case( excit_type%double_L_to_R )
        case( excit_type%double_R_to_L )
        case( excit_type%fullstop_lowering )
        case( excit_type%fullstop_raising )
        case( excit_type%fullstop_L_to_R )
        case( excit_type%fullstop_R_to_L )
        case( excit_type%fullstart_lowering)
        case( excit_type%fullstart_raising)
        case( excit_type%fullstart_L_to_R)
        case( excit_type%fullstart_R_to_L)
        case( excit_type%fullstart_stop_alike)
        case( excit_type%fullstart_stop_mixed)
        case default
            print *, "incorrect typ: ", excit_names(typ)
            call stop_all(this_routine, "see above")
        end select
#endif

        ASSERT(all(inds > 0) .and. all(inds <= nSpatOrbs))

        excit_info_int = 0_int64

        call encode_excit_info_type(excit_info_int, typ)
        call encode_excit_info_indices(excit_info_int, inds)

    end function encode_excit_info_vec

    function encode_excit_info_scalar(typ, a, i, b, j) result(excit_info_int)
        ! function to encode the minimal information of an excit-info
        ! object into a single 64bit integer. used in the PCHB excitation
        ! generation.
        debug_function_name("encode_excit_info_scalar")
        integer, intent(in) :: typ, a, i, b, j
        integer(int64) :: excit_info_int

#ifdef DEBUG_
        select case(typ)
        case( excit_type%single_overlap_L_to_R)
        case( excit_type%single_overlap_R_to_L )
        case( excit_type%double_lowering )
        case( excit_type%double_raising )
        case( excit_type%double_L_to_R_to_L)
        case( excit_type%double_R_to_L_to_R )
        case( excit_type%double_L_to_R )
        case( excit_type%double_R_to_L )
        case( excit_type%fullstop_lowering )
        case( excit_type%fullstop_raising )
        case( excit_type%fullstop_L_to_R )
        case( excit_type%fullstop_R_to_L )
        case( excit_type%fullstart_lowering)
        case( excit_type%fullstart_raising)
        case( excit_type%fullstart_L_to_R)
        case( excit_type%fullstart_R_to_L)
        case( excit_type%fullstart_stop_alike)
        case( excit_type%fullstart_stop_mixed)
        case default
            print *, "incorrect typ: ", excit_names(typ)
            call stop_all(this_routine, "see above")
        end select
#endif

        ASSERT(a > 0 .and. a <= nSpatOrbs)
        ASSERT(i > 0 .and. i <= nSpatOrbs)
        ASSERT(b > 0 .and. b <= nSpatOrbs)
        ASSERT(j > 0 .and. j <= nSpatOrbs)

        excit_info_int = encode_excit_info_vec(typ, [a,i,b,j])

    end function encode_excit_info_scalar

    subroutine encode_excit_info_type(excit_info_int, typ)
        debug_function_name("encode_excit_info_type")
        integer(int64), intent(inout) :: excit_info_int
        integer, intent(in) :: typ

#ifdef DEBUG_
        select case(typ)
        case( excit_type%single_overlap_L_to_R)
        case( excit_type%single_overlap_R_to_L )
        case( excit_type%double_lowering )
        case( excit_type%double_raising )
        case( excit_type%double_L_to_R_to_L)
        case( excit_type%double_R_to_L_to_R )
        case( excit_type%double_L_to_R )
        case( excit_type%double_R_to_L )
        case( excit_type%fullstop_lowering )
        case( excit_type%fullstop_raising )
        case( excit_type%fullstop_L_to_R )
        case( excit_type%fullstop_R_to_L )
        case( excit_type%fullstart_lowering)
        case( excit_type%fullstart_raising)
        case( excit_type%fullstart_L_to_R)
        case( excit_type%fullstart_R_to_L)
        case( excit_type%fullstart_stop_alike)
        case( excit_type%fullstart_stop_mixed)
        case default
            print *, "incorrect typ: ", excit_names(typ)
            call stop_all(this_routine, "see above")
        end select
#endif

        ! the convention is to store the excit-type 'first' at the LSB
        ! so this should be fine:
        call mvbits(int(typ, int64), 0, n_excit_type_bits, excit_info_int, 0)

    end subroutine encode_excit_info_type

    pure function extract_excit_info_type(excit_info_int) result(typ)
        integer(int64), intent(in) :: excit_info_int
        integer :: typ

        typ = int(ibits(excit_info_int, 0, n_excit_type_bits))

    end function extract_excit_info_type

    subroutine encode_excit_info_indices_vec(excit_info_int, inds)
        debug_function_name("encode_excit_info_indices_vec")
        integer(int64), intent(inout) :: excit_info_int
        integer, intent(in) :: inds(4)

        integer :: i

        ! i already used 5 bits for the excit-type: so 64-5 = 59 remaining
        ! 59/4 = 14.. so theoretically 14 bits remaining for each SPATIAL
        ! orbital.. 2^14 > 16k orbitals... not necessary.. maybe store more
        ! info. for now i use 8 bits: 256 SPATIAL orbital max.. that
        ! should be enough..
        ASSERT(all(inds > 0) .and. all(inds <= nSpatOrbs))
        ASSERT(nSpatOrbs <= 256)

        do i = 1, 4
            call mvbits(int(inds(i), int64), 0, n_excit_index_bits, &
                excit_info_int, n_excit_type_bits + (i - 1) * n_excit_index_bits)
        end do

    end subroutine encode_excit_info_indices_vec

    subroutine encode_excit_info_indices_scalar(excit_info_int, a, i, b, j)
        integer(int64), intent(inout) :: excit_info_int
        integer, intent(in) :: a, i, b, j

        call encode_excit_info_indices_vec(excit_info_int, [a,i,b,j])

    end subroutine encode_excit_info_indices_scalar

    subroutine extract_excit_info_indices_vec(excit_info_int, inds)
        integer(int64), intent(in) :: excit_info_int
        integer, intent(out) :: inds(4)

        integer :: i

        inds = 0
        do i = 1, 4
            inds(i) = extract_excit_info_index(excit_info_int, i)
        end do

    end subroutine extract_excit_info_indices_vec

    subroutine extract_excit_info_indices_scalar(excit_info_int, a, i, b, j)
        integer(int64), intent(in) :: excit_info_int
        integer, intent(out) :: a, i, b, j

        a = extract_excit_info_index(excit_info_int, 1)
        i = extract_excit_info_index(excit_info_int, 2)
        b = extract_excit_info_index(excit_info_int, 3)
        j = extract_excit_info_index(excit_info_int, 4)

    end subroutine extract_excit_info_indices_scalar

    function extract_excit_info_index(excit_info_int, pos) result(ind)
        debug_function_name("extract_excit_info_index")
        integer(int64), intent(in) :: excit_info_int
        integer, intent(in) :: pos
        integer :: ind

        ASSERT(pos > 0 .and. pos <= 4)

        ind = int(ibits(excit_info_int, n_excit_type_bits &
            + (pos - 1) * n_excit_index_bits, n_excit_index_bits))

    end function extract_excit_info_index

    subroutine extract_excit_info_scalar(excit_info_int, typ, a, i, b, j)
        integer(int64), intent(in) :: excit_info_int
        integer, intent(out) :: typ, a, i, b, j

        typ = extract_excit_info_type(excit_info_int)

        call extract_excit_info_indices_scalar(excit_info_int, a, i, b, j)

    end subroutine extract_excit_info_scalar

    subroutine extract_excit_info_vector(excit_info_int, typ, inds)
        integer(int64), intent(in) :: excit_info_int
        integer, intent(out) :: typ, inds(4)

        typ = extract_excit_info_type(excit_info_int)

        call extract_excit_info_indices_vec(excit_info_int, inds)

    end subroutine extract_excit_info_vector

    subroutine extract_excit_info_obj(excit_info_int, excitInfo)
        integer(int64), intent(in) :: excit_info_int
        type(ExcitationInformation_t), intent(out) :: excitInfo

        integer :: typ, a, i, b, j

        ! this function needs to do additional processing of the
        ! minimal info stored in excit_info_int

        typ = extract_excit_info_type(excit_info_int)

        call extract_excit_info_indices_scalar(excit_info_int, a, i, b, j)

        excitInfo = calc_remaining_excit_info(typ, a, i, b, j)

    end subroutine extract_excit_info_obj

    function calc_remaining_excit_info(typ, a, i, b, j) result(excitInfo)
        debug_function_name("calc_remaining_excit_info")
        integer, intent(in) :: typ, a, i, b, j
        type(ExcitationInformation_t) :: excitInfo

        integer :: gen1, gen2, currentGen, firstGen, lastGen, fullstart, &
                   secondStart, firstEnd, fullEnd, weight, excitLvl, overlap
        real(dp) :: order, order1

#ifdef DEBUG_
        select case(typ)
        case( excit_type%single_overlap_L_to_R)
        case( excit_type%single_overlap_R_to_L )
        case( excit_type%double_lowering )
        case( excit_type%double_raising )
        case( excit_type%double_L_to_R_to_L)
        case( excit_type%double_R_to_L_to_R )
        case( excit_type%double_L_to_R )
        case( excit_type%double_R_to_L )
        case( excit_type%fullstop_lowering )
        case( excit_type%fullstop_raising )
        case( excit_type%fullstop_L_to_R )
        case( excit_type%fullstop_R_to_L )
        case( excit_type%fullstart_lowering)
        case( excit_type%fullstart_raising)
        case( excit_type%fullstart_L_to_R)
        case( excit_type%fullstart_R_to_L)
        case( excit_type%fullstart_stop_alike)
        case( excit_type%fullstart_stop_mixed)
        case default
            print *, "incorrect typ: ", excit_names(typ)
            call stop_all(this_routine, "see above")
        end select
#endif

        ASSERT(a > 0 .and. a <= nSpatOrbs)
        ASSERT(i > 0 .and. i <= nSpatOrbs)
        ASSERT(b > 0 .and. b <= nSpatOrbs)
        ASSERT(j > 0 .and. j <= nSpatOrbs)

        ! do the necessary recomputation of excitInfo entries
        ! in Debug mode also check if the indices are correct!
        ! and for now assume that this function gets only called in the
        ! pchb list creation, which quite restricts the indices.
        ! this might need to be changed in the future, but for now it is
        ! good to also ensure the PCHB lists are created correctly!

        ! set up some defaults which are only changed if necessary:

        weight = 0
        ! fuck this excitLvl info is so stupid... i need to change that :(
        excitLvl = 2
        order = 1.0_dp
        order1 = 1.0_dp

        select case (typ)

        case ( excit_type%single_overlap_L_to_R )

            ASSERT(a == b)
            ASSERT(i < j)
            ASSERT(a > i .and. a < j)

            gen1        = gen_type%L
            gen2        = gen_type%R
            currentGen  = gen_type%L
            firstGen    = gen_type%L
            lastGen     = gen_type%R
            fullstart   = i
            secondStart = a
            firstEnd    = a
            fullEnd     = j
            overlap     = 1



        case ( excit_type%single_overlap_R_to_L )

            ASSERT(i == j)
            ASSERT(a < b)
            ASSERT( i > a .and. i < b)

            gen1        = gen_type%R
            gen2        = gen_type%L
            currentGen  = gen_type%R
            firstGen    = gen_type%R
            lastGen     = gen_type%L
            fullstart   = a
            secondStart = i
            firstEnd    = i
            fullEnd     = b
            overlap     = 1

        case ( excit_type%double_lowering )

            ASSERT(i < j)
            ASSERT(j < a)
            ASSERT(a < b)

            gen1        = gen_type%L
            gen2        = gen_type%L
            currentGen  = gen_type%L
            firstGen    = gen_type%L
            lastGen     = gen_type%L
            fullstart   = i
            secondStart = j
            firstEnd    = a
            fullEnd     = b
            overlap     = firstEnd - secondStart + 1


        case ( excit_type%double_raising )

            ! in this assignment i switch to E_{bj}E_{ai} assignement to
            ! ensure 'correct' sign convention of the x1 coupling coeffs!
            ASSERT(b < a)
            ASSERT(a < j)
            ASSERT(j < i)

            gen1        = gen_type%R
            gen2        = gen_type%R
            currentGen  = gen_type%R
            firstGen    = gen_type%R
            lastGen     = gen_type%R
            fullstart   = b
            secondStart = a
            firstEnd    = j
            fullEnd     = i

            overlap     = firstEnd - secondStart + 1

        case ( excit_type%double_L_to_R_to_L )
            ! here we switch to E_{aj}E_{bi}

            ASSERT( j < a )
            ASSERT( a < i )
            ASSERT( i < b )

            gen1        = gen_type%R
            gen2        = gen_type%L
            currentGen  = gen_type%L
            firstGen    = gen_type%L
            lastGen     = gen_type%L
            fullstart   = j
            secondStart = a
            firstEnd    = i
            fullEnd     = b

            overlap     = firstEnd - secondStart + 1


        case ( excit_type%double_R_to_L_to_R )

            ! switch to E_{aj}E_{bi}
            ASSERT( a < j )
            ASSERT( j < b )
            ASSERT( b < i )

            gen1        = gen_type%R
            gen2        = gen_type%L
            currentGen  = gen_type%R
            firstGen    = gen_type%R
            lastGen     = gen_type%R
            fullstart   = a
            secondStart = j
            firstEnd    = b
            fullEnd     = i

            overlap = firstEnd - secondStart + 1

        case ( excit_type%double_L_to_R )

            ! switch to E_{aj}E_{bi}
            ASSERT( j < a )
            ASSERT( a < b )
            ASSERT( b < i )

            gen1        = gen_type%R
            gen2        = gen_type%L
            currentGen  = gen_type%L
            firstGen    = gen_type%L
            lastGen     = gen_type%R
            fullstart   = j
            secondStart = a
            firstEnd    = b
            fullEnd     = i

            overlap = firstEnd - secondStart + 1


        case ( excit_type%double_R_to_L )

            ! switch to E_{aj}E_{bi}
            ASSERT( a < j )
            ASSERT( j < i )
            ASSERT( i < b )

            gen1        = gen_type%R
            gen2        = gen_type%L
            currentGen  = gen_type%R
            firstGen    = gen_type%R
            lastGen     = gen_type%L
            fullstart   = a
            secondStart = j
            firstEnd    = i
            fullEnd     = b

            overlap = firstEnd - secondStart + 1


        case ( excit_type%fullstop_lowering )

            ASSERT( i < j )
            ASSERT( j < a )
            ASSERT( a == b )

            gen1        = gen_type%L
            gen2        = gen_type%L
            currentGen  = gen_type%L
            firstGen    = gen_type%L
            lastGen     = gen_type%L
            fullstart   = i
            secondStart = j
            firstEnd    = a
            fullEnd     = a

            overlap = firstEnd - secondStart + 1


        case ( excit_type%fullstop_raising )

            ASSERT(a < b )
            ASSERT(b < i )
            ASSERT(i == j)

            gen1        = gen_type%R
            gen2        = gen_type%R
            currentGen  = gen_type%R
            firstGen    = gen_type%R
            lastGen     = gen_type%R
            fullstart   = a
            secondStart = b
            firstEnd    = i
            fullEnd     = i

            overlap = firstEnd - secondStart + 1


        case ( excit_type%fullstop_R_to_L )

            ! switch to E_{aj}E_{bi}
            ASSERT(a < j )
            ASSERT(j < b )
            ASSERT(b == i )

            gen1        = gen_type%R
            gen2        = gen_type%L
            currentGen  = gen_type%R
            firstGen    = gen_type%R
            lastGen     = gen_type%L
            fullstart   = a
            secondStart = j
            firstEnd    = b
            fullEnd     = b

            overlap = firstEnd - secondStart + 1


        case ( excit_type%fullstop_L_to_R )

            ! switch to E_{aj}E_{bi}
            ASSERT( j < a )
            ASSERT( a < b )
            ASSERT( b == i )

            gen1        = gen_type%R
            gen2        = gen_type%L
            currentGen  = gen_type%L
            firstGen    = gen_type%L
            lastGen     = gen_type%R
            fullstart   = j
            secondStart = a
            firstEnd    = b
            fullEnd     = b

            overlap = firstEnd - secondStart + 1


        case ( excit_type%fullstart_lowering )

            ASSERT( i == j)
            ASSERT( j < a )
            ASSERT( a < b )

            gen1        = gen_type%L
            gen2        = gen_type%L
            currentGen  = gen_type%L
            firstGen    = gen_type%L
            lastGen     = gen_type%L
            fullstart   = i
            secondStart = i
            firstEnd    = a
            fullEnd     = b

            overlap = firstEnd - secondStart + 1


        case ( excit_type%fullstart_raising )

            ASSERT( a == b )
            ASSERT( b < i )
            ASSERT( i < j )

            gen1        = gen_type%R
            gen2        = gen_type%R
            currentGen  = gen_type%R
            firstGen    = gen_type%R
            lastGen     = gen_type%R
            fullstart   = a
            secondStart = a
            firstEnd    = i
            fullEnd     = j

            overlap = firstEnd - secondStart + 1


        case ( excit_type%fullstart_L_to_R )

            ! switch to E_{aj} E_{bi}
            ASSERT( a == j )
            ASSERT( j < b )
            ASSERT( b < i )

            gen1        = gen_type%R
            gen2        = gen_type%L
            currentGen  = gen_type%R
            firstGen    = gen_type%L
            lastGen     = gen_type%R
            fullstart   = a
            secondStart = a
            firstEnd    = b
            fullEnd     = i

            overlap = firstEnd - secondStart + 1


        case ( excit_type%fullstart_R_to_L )

            ! switch to E_{aj}E_{bi}
            ASSERT( a == j )
            ASSERT( j < i )
            ASSERT( i < b )

            gen1        = gen_type%R
            gen2        = gen_type%L
            currentGen  = gen_type%L
            firstGen    = gen_type%R
            lastGen     = gen_type%L
            fullstart   = a
            secondStart = a
            firstEnd    = i
            fullEnd     = b

            overlap = firstEnd - secondStart + 1


        case ( excit_type%fullstart_stop_alike )

            ! here i need if statement..

            ASSERT( i == j )
            ASSERT( a == b )
            ASSERT( a /= i )

            if (a < i) then
                gen1 = gen_type%R
            else if (a > i) then
                gen1 = gen_type%L
            end if

            gen2        = gen1
            currentGen  = gen1
            firstGen    = gen1
            lastGen     = gen1
            fullstart   = min(i,a)
            secondStart = min(i,a)
            firstEnd    = max(i,a)
            fullEnd     = max(i,a)

            overlap = firstEnd - secondStart + 1

        case ( excit_type%fullstart_stop_mixed )

            ! switch to E_{aj}E_{bi}
            ASSERT(a == j)
            ASSERT(b == i)
            ASSERT(a /= b)

            gen1        = gen_type%R
            gen2        = gen_type%L
            currentGen  = gen_type%R
            firstGen    = gen_type%R
            lastGen     = gen_type%R
            fullstart   = min(a,i)
            secondStart = min(a,i)
            firstEnd    = max(a,i)
            fullEnd     = max(a,i)

            overlap = firstEnd - secondStart + 1

        end select

        excitInfo = assign_excitInfo_values_exact(typ, gen1, gen2, currentGen, &
            firstGen, lastGen, a, i, b, j, fullstart, secondStart, firstEnd, &
            fullEnd, weight, excitLvl, order, order1, overlap)

    end function calc_remaining_excit_info

    function encode_excit_info_obj(excitInfo) result(excit_info_int)
        type(ExcitationInformation_t), intent(in) :: excitInfo
        integer(int64) :: excit_info_int

        ! function to encode directly from ExcitationInformation_t type to
        ! a 64-bit integer
        excit_info_int = 0_int64

        excit_info_int =  encode_excit_info_vec(excitInfo%typ, &
            [excitInfo%i, excitInfo%j, excitInfo%k, excitInfo%l])

    end function encode_excit_info_obj

end module guga_bitRepOps
