#include "macros.h"

! create a new module which covers all the real-space hubbard stuff in a more 
! clean and more efficient way. 
! I also want to implement triangular and kagome lattice types, so it is 
! good to get this better sorted. 
! also i definetly have to improve on the hubbard excitation generator
! i also probably should decouple all the UEG k-space hubbard at a certain 
! point. because it makes things super messy and also more inefficient how it 
! is done right now..
! and for the beginning i should probably stick to the already working 
! lattices in NECI
! and i also want to make it compatible to restart old runs with this new 
! implementation

module real_space_hubbard

    use SystemData, only: t_new_real_space_hubbard, lattice_type, length_x, &
                          length_y, length_z, uhub, nbasis, bhub, t_open_bc_x, &
                          t_open_bc_y, t_open_bc_z, G1, ecore, nel, nOccAlpha, nOccBeta, &
                          t_trans_corr, trans_corr_param, t_trans_corr_2body, & 
                          trans_corr_param_2body, tHPHF, t_trans_corr_new
    use lattice_mod, only: lattice, determine_optimal_time_step, lat, &
                    get_helement_lattice, get_helement_lattice_ex_mat, & 
                    get_helement_lattice_general
    use constants, only: dp, EPS, n_int, bits_n_int
    use procedure_pointers, only: get_umat_el, generate_excitation
    use OneEInts, only: tmat2d, GetTMatEl
    use fcimcdata, only: pSingles, pDoubles, tsearchtau, tsearchtauoption
    use CalcData, only: t_hist_tau_search, t_hist_tau_search_option, tau
!     use umatcache, only: gtid
    use dsfmt_interface, only: genrand_real2_dsfmt
    use DetBitOps, only: FindBitExcitLevel, EncodeBitDet
    use bit_rep_data, only: NIfTot
    use util_mod, only: binary_search_first_ge, choose
    use bit_reps, only: decode_bit_det
    use sort_mod, only: sort
!     use fcimc_helper, only: update_run_reference
!     use Determinants, only: write_det

    implicit none 

    real(dp) :: lat_tau_factor = 0.5_dp

    ! create a flag which indicate to start in a neel state 
    logical :: t_start_neel_state = .false. 

    interface get_helement_rs_hub
        module procedure get_helement_rs_hub_ex_mat
        module procedure get_helement_rs_hub_general
    end interface get_helement_rs_hub

contains 

    ! some brainstorming: 
    ! i want to change the input so that one specifies the lattice type 
    ! in the System input by 
    ! system hubbard [lattice-type] 
    subroutine init_real_space_hubbard() 
        use SystemData, only: tExch, thub, treal
        ! routine, which does all of the necessary initialisation
        character(*), parameter :: this_routine = "init_real_space_hubbard"
        ! especially do some stop_all here so no wrong input is used 
        real(dp) :: tau_opt
        integer :: neel_state_ni(nel)
        integer(n_int) :: ilut_neel(0:NIfTot)

        print *, "using new real-space hubbard implementation: "


        ! i do not need exchange integrals in the real-space hubbard model
        tExch = .false.
        ! after the whole setup i can set thub to false or? 
        thub = .false.
        ! and treal i can also set to false or? 
        treal = .false.
        ! first assert all the right input! 
        call check_real_space_hubbard_input() 

        ! which stuff do i need to initialize here? 
        get_umat_el => get_umat_el_hub

        ! also use the new lattice matrix elements

        ! i have to check if the lattice should be constructed from an fcidump 
        ! or created internally.. 
        if (trim(adjustl(lattice_type)) == 'read') then 
            ! then i have to construct tmat first 
            call init_tmat() 
            ! and then construct the lattice 
            lat => lattice(lattice_type, length_x, length_y, length_z, .not. t_open_bc_x, &
                .not. t_open_bc_y, .not. t_open_bc_z)
        else 
            ! otherwise i have to do it the other way around 
            lat => lattice(lattice_type, length_x, length_y, length_z, .not. t_open_bc_x, &
                .not. t_open_bc_y, .not. t_open_bc_z)
     
            ! if nbaiss was not yet provided:
            if (nbasis <= 0) then 
                nbasis = 2 * lat%get_nsites() 
            end if

            call init_tmat(lat)

        end if

        ! i guess i have to setup G1 also.. argh.. i hate this! 
        allocate(G1(nbasis)) 
        G1(1:nbasis-1:2)%ms = -1
        G1(2:nbasis:2)%ms = 1

        ! Ecore should default to 0, but be sure anyway! 
        ecore = 0.0_dp

        ! and i have to point to the new hubbard excitation generator
        pSingles = 1.0_dp 
        pDoubles = 0.0_dp

        ! and i have to calculate the optimal time-step for the hubbard models. 
        ! where i need the connectivity of the lattice i guess? 
        if (.not. tHPHF) then
            generate_excitation => gen_excit_rs_hubbard
        end if
        
        ! i have to calculate the optimal time-step
        ! and maybe i have to be a bit more safe here and not be too near to 
        ! the optimal time-step
        tau_opt = determine_optimal_time_step()
        if (tau < EPS) then 
            print *, "setting time-step to optimally determined time-step: ", tau_opt
            print *, "times: ", lat_tau_factor
            tau = lat_tau_factor * tau_opt

        else 
            print *, "optimal time-step would be: ", tau_opt
            print *, "but tau specified in input!"
        end if

        ! and i have to turn off the time-step search for the hubbard 
        tsearchtau = .false.
        ! set tsearchtauoption to true to use the death-tau search option
        tsearchtauoption = .true.

        t_hist_tau_search = .false. 
        t_hist_tau_search_option = .false. 
        
        if (t_start_neel_state) then 
!             neel_state_ni = create_neel_state(ilut_neel)

            print *, "starting from the Neel state: " 
            if (nel > nbasis/2) then 
                call stop_all(this_routine, &
                    "more than half-filling! does neel state make sense?")
            end if
!             call write_det(6, neel_state_ni, .true.)
!             call changerefdet(neel_state_ni)
!             call update_run_reference(ilut_neel, 1)

        end if
        ! do not set that here, due to circular dependencies
!         max_death_cpt = 0.0_dp

        call init_get_helement_hubbard()

    end subroutine init_real_space_hubbard

    subroutine init_get_helement_hubbard
        get_helement_lattice_ex_mat => get_helement_rs_hub_ex_mat
        get_helement_lattice_general => get_helement_rs_hub_general
        call init_tmat(lat)
    end subroutine init_get_helement_hubbard

    subroutine check_real_space_hubbard_input() 
        use SystemData, only: tCSF, tReltvy, tUEG, tUEG2, tHub, & 
                              tKPntSym, tLatticeGens, tUEGNewGenerator, &
                tGenHelWeighted, tGen_4ind_weighted, tGen_4ind_reverse, &
                tUEGNewGenerator, tGen_4ind_part_exact, tGen_4ind_lin_exact, &
                tGen_4ind_2, tGen_4ind_2_symmetric, tGen_4ind_unbound, tStoreSpinOrbs, &
                tReal
!         use umatcache, only : tTransGTid
        use OneEInts, only: tcpmdsymtmat, tOneelecdiag

        character(*), parameter :: this_routine = "check_real_space_hubbard_input"
        ! do all the input checking here, so no wrong input is used!

        if (tCSF)             call stop_all(this_routine, "tCSF set to true!")
        if (tReltvy)          call stop_all(this_routine, "tReltvy set to true!")

        ! what else.. 
        if (tUEG)             call stop_all(this_routine, "tUEG set to true!")
        if (tUEG2)            call stop_all(this_routine, "tUEG2 set to true!")
        if (tHub)             call stop_all(this_routine, "tHub set to true!")
        if (tReal)            call stop_all(this_routine, "tReal set to true!")
        if (tKPntSym)         call stop_all(this_routine, "tKPntSym set to true!")
        if (tLatticeGens)     call stop_all(this_routine, "tLatticeGens set to true!")
        if (tUEGNewGenerator) call stop_all(this_routine, "tUEGNewGenerator set to true!")
        if (tGenHelWeighted)  call stop_all(this_routine, "tGenHelWeighted") 
        if (tGen_4ind_weighted) call stop_all(this_routine, "tGen_4ind_weighted") 
        if (tGen_4ind_reverse) call stop_all(this_routine, "tGen_4ind_reverse") 
        if (tGen_4ind_part_exact) call stop_all(this_routine, "tGen_4ind_part_exact") 
        if (tGen_4ind_2)        call stop_all(this_routine, "tGen_4ind_2") 
        if (tGen_4ind_2_symmetric) call stop_all(this_routine, "tGen_4ind_2_symmetric") 
        if (tGen_4ind_unbound)      call stop_all(this_routine, "tGen_4ind_unbound")
        if (tStoreSpinOrbs)     call stop_all(this_routine, "tStoreSpinOrbs")
!         if (tTransGTid)         call stop_all(this_routine, "tTransGTid")
        if (tcpmdsymtmat)        call stop_all(this_routine, "tcpmdsymmat")
        if (tOneelecdiag)       call stop_all(this_routine, "tOneelecdiag")
            
    end subroutine check_real_space_hubbard_input

    ! i could also create all determinants, not only the open-shells.. or? 
    function create_all_dets(n_orbs, n_alpha, n_beta) result(all_dets)
        integer, intent(in) :: n_orbs, n_alpha, n_beta
        integer(n_int), allocatable :: all_dets(:) 
        character(*), parameter :: this_routine = "create_all_dets" 

        integer, allocatable :: n_dets(:)
        integer(n_int), allocatable :: alpha_basis(:), beta_basis(:)
        integer :: i, n_states, n_total, j

        n_dets = calc_n_double(n_orbs, n_alpha, n_beta) 

        ! cannot deal with more than 32 orbitals yet.. since i am lazy! 
        ASSERT(n_orbs <= 32) 

        ! do it in a really simple way for now.. 
        ! just run over all the alpha and beta spin-bases and combine them 
        ! to the full-basis 

        alpha_basis = create_one_spin_basis(n_orbs, n_alpha) 
        beta_basis = create_one_spin_basis(n_orbs, n_beta) 

        ! keep an count on all the states
        n_states = 0

        ! and allocate the array, which keeps all the states 
        n_total = sum(n_dets) 

        allocate(all_dets(n_total))
        all_dets = 0_n_int 

        ! check if everything went correctly: 
        ASSERT(size(alpha_basis) * size(beta_basis) == n_total)

        do i = 1, size(alpha_basis) 
            do j = 1, size(beta_basis)
                n_states = n_states + 1
                all_dets(n_states) = general_product(alpha_basis(i),beta_basis(j),n_orbs)
            end do
        end do
        
        call sort(all_dets)

        ASSERT(n_states == n_total) 

    end function create_all_dets 

    subroutine create_hilbert_space_realspace(n_orbs, n_alpha, n_beta, & 
            n_states, state_list_ni, state_list_ilut) 
        integer, intent(in) :: n_orbs, n_alpha, n_beta
        integer, intent(out) :: n_states 
        integer, intent(out), allocatable :: state_list_ni(:,:)
        integer(n_int), intent(out), allocatable :: state_list_ilut(:,:) 

        integer(n_int), allocatable :: all_dets(:) 
        integer(n_int) :: temp_ilut(0:niftot)
        integer :: nJ(nel), i

        all_dets = create_all_dets(n_orbs, n_alpha, n_beta) 

        n_states = size(all_dets)

        allocate(state_list_ni(nel,n_states))
        allocate(state_list_ilut(0:niftot,n_states))

        do i = 1, size(all_dets) 
            temp_ilut = all_dets(i) 
            call decode_bit_det(nJ, temp_ilut)

            state_list_ni(:,i) = nJ
            state_list_ilut(:,i) = temp_ilut

        end do

    end subroutine create_hilbert_space_realspace

    function create_all_open_shell_dets(n_orbs, n_alpha, n_beta) result(open_shells)
        integer, intent(in) :: n_orbs, n_alpha, n_beta
        ! Alis idea for the use of the ADI option for the real-space hubbard 
        ! is to let all fully open shell dets, and excitations thereof 
        ! (dets with exactly one double occupancy in the case of half-filling) 
        ! be initiators 
        character(*), parameter :: this_routine = "create_all_open_shell_dets" 

        integer, allocatable :: n_dets(:)
        integer(n_int), allocatable :: open_shells(:), max_basis(:)
        integer :: n_max, n_min

        n_dets = calc_n_double(n_orbs, n_alpha, n_beta)

        ! only implement that for less than 32 orbitals, since i do not want to 
        ! deal with multiple integers to store the basis and also the size of 
        ! the hilbert space would be too big anyway.. 
        ASSERT(n_orbs <= 32)

        ! how should i create all of those purely open-shell dets? 
        ! in the case of half-filling it is easiest.. 
        
        ! i could first distribute the alpha spins and then the beta spins 
        ! afterwards.. 
        ! i should work with the bit-representation, then i get to know the 
        ! fortran intrinsic bit operations again.. 

        n_max = max(n_alpha, n_beta)
        n_min = min(n_alpha, n_beta)
        
        ! first create the spin-basis with more electrons: 
        max_basis = create_one_spin_basis(n_orbs, n_max) 

        open_shells = combine_spin_basis(n_orbs, n_max, n_min, n_dets(1),  max_basis, .true.)

    end function create_all_open_shell_dets

    function combine_spin_basis(n_orbs, n_first, n_second, n_total, first_basis, & 
            t_sort, n_doubles) result(spin_basis) 
        ! function to combine a already given spin-basis with the second one
        ! with the additional option to specify the number of doubly occupied 
        ! determinants.. 
        ! n_orbs   ...  is the number of spatial orbitals
        ! n_first  ...  number of spins, with which first_basis was created 
        ! n_second ...  number of spins with which basis is combined now 
        ! n_total  ...  total size of basis for given number of doubly occ. sites!
        ! first_basis . already contructed one-spin basis
        ! t_sort   ...  should the basis be sorted. (after combination i am not  sure it will be sorted..
        ! n_doubles ..  (optional:) number of doubly occupied sites(has to fit with n_total!)
        integer, intent(in) :: n_orbs, n_first, n_second, n_total
        integer(n_int), intent(in) :: first_basis(:) 
        logical, intent(in) :: t_sort
        integer, intent(in), optional :: n_doubles 
        integer(n_int) :: spin_basis(n_total) 
        character(*), parameter :: this_routine = "combine_spin_basis" 
#ifdef __DEBUG
        integer, allocatable :: n_dets(:)
#endif 
        integer :: n_doub, i, j, n_remain, n
        integer(n_int), allocatable :: second_basis(:)
    
        ! the half-filled case is the most easy one.. should we treat it 
        ! special? maybe.. 

        ASSERT(n_first >= n_second) 
#ifdef __DEBUG 
        ! be sure that the provided n_total fits with the n_doubles if 
        ! provided 
        n_dets = calc_n_double(n_orbs, n_first, n_second)
        if (present(n_doubles)) then 
            ASSERT(n_dets(n_doubles+1) == n_total)
        else 
            ! n_doubles = 0 is the default
            ASSERT(n_dets(1) == n_total)
        end if
#endif

        ! default value of n_doubles = 0 
        if (present(n_doubles)) then 
            n_doub = n_doubles
        else 
            n_doub = 0 
        end if

        ! do the n_double = 0 and half-filled case first, thats the 
        ! easiest and only necessary for now.. 
        select case (n_doub) 
        case (0) 

            if (n_first + n_second == n_orbs) then 
                ! half-filled case: easiest 
                ! since ms is not important, just convert the already 
                ! calculated spin-basis to a spin-orbital basis: 
                ! 0011 -> 00 00 01 01 eg. 
                do i = 1, n_total 
                    ! do have it ordered -> first set the beta spins on the left
                    spin_basis(i) = set_alpha_beta_spins(first_basis(i), n_orbs,.false.) 
                    ! and combine it with the alpha spins..
                    spin_basis(i) = ieor(spin_basis(i),set_alpha_beta_spins(not(first_basis(i)), n_orbs,.true.))
                end do
            else 
                ! we have to distribute the n_second remaining spins 
                ! across the n_orbs - n_first orbitals, so there are: 
                n_remain = int(choose(n_orbs - n_first, n_second))
                ! states per first_basis states 
                ASSERT(n_remain * size(first_basis) == n_total) 

                second_basis = create_one_spin_basis(n_orbs - n_first, n_second)

                ! i want to do some sort of tensor product here:
                ! |0011> x |01> = |00 10 01 01>
                ! |0011> x |10> = |10 00 01 01>  
                n = 1 
                do i = 1, size(first_basis) 
                    do j = 1, size(second_basis)
                        spin_basis(n) = open_shell_product(first_basis(i),second_basis(j),n_orbs)
                        n = n + 1
                    end do
                end do

            end if

        case default 

            ! i should not need a select case here.. this gets to 
            ! complicated.. 
            ! it is similar to non-half filled actually.. 
            if (n_first + n_second == n_orbs) then 
                ! then we have the half-filled case.. 
            else

            end if

            ! not yet implemented 
            call stop_all(this_routine, "not yet implemented!")

        end select

        ! do i want to sort it?? 
        if (t_sort) then 
            call sort(spin_basis) 
        end if


    end function combine_spin_basis

    function general_product(alpha, beta, n_orbs) result(basis) 
        ! this is a "general" tensor product of two spin basisfunctions.
        ! the ones in the alpha integer must be set at position 2*i 
        ! and the beta at 2*j-1 
        integer(n_int), intent(in) :: alpha, beta 
        integer, intent(in) :: n_orbs 
        integer(n_int) :: basis 

        ! it should be as simple: 
        basis = set_alpha_beta_spins(alpha, n_orbs, .false.) 
        basis = xor(basis, set_alpha_beta_spins(beta, n_orbs, .true.))

    end function general_product

    function open_shell_product(alpha, beta, n_orbs) result(basis)
        ! compute the tensor product of alpha and beta spin-bases to give 
        ! fully open-shell determinants, the first input is to be defined as 
        ! the alpha spin part
        ! and the beta input then is only as big as n_orbs - n_alpha and only 
        ! tells us in which empty orbitals of the alpha orbitals beta spins 
        ! should be set 
        ! sorting is not so easily ensured here anymore.. 
        integer(n_int), intent(in) :: alpha, beta
        integer, intent(in) :: n_orbs
        integer(n_int) :: basis

        integer, allocatable :: nZeros(:), nOnes(:)
        integer(n_int) :: mask_zeros
        integer :: i
        ! the setting of the alpha spins is the same or? 
        basis = set_alpha_beta_spins(alpha, n_orbs, .false.)

        ! i need the zeros, but just in the n_orbs range
        mask_zeros = iand(not(alpha), maskr(n_orbs))

        allocate(nZeros(popcnt(mask_zeros)))

        call decode_bit_det(nZeros, [mask_zeros])

        ! now i need the ones from beta
        allocate(nOnes(popcnt(beta)))
        call decode_bit_det(nOnes, [beta])

        ! and i have to set all beta-spins indicated by nZeros(nOnes) 
        ! so we need beta spin! 
        nZeros = 2*nZeros - 1
        do i = 1, size(nOnes) 
            basis = ibset(basis, nZeros(nOnes(i))-1)
        end do

    end function open_shell_product

    function set_alpha_beta_spins(beta_mask, n_orbs, t_beta) result(beta_spins) 
        ! a function which converts a flag of spatial beta-spins to a 
        ! spin orbita basis, eg.: 0011 -> 00 00 01 01
        integer(n_int), intent(in) :: beta_mask
        integer, intent(in) :: n_orbs 
        logical, intent(in) :: t_beta
        integer(n_int) :: beta_spins

        integer, allocatable :: nOnes(:) 
        integer :: i 

        allocate(nOnes(popcnt(iand(beta_mask, maskr(n_orbs))))) 

        call decode_bit_det(nOnes, [beta_mask])

        if (t_beta) then 
            ! then we want to set beta spins: 
            nOnes = 2*nOnes - 1
        else 
            ! otherwise we want to set alpha spins
            nOnes = 2*nOnes
        end if
                
        beta_spins = 0_n_int
        do i = 1, size(nOnes) 
            beta_spins = ibset(beta_spins, nOnes(i)-1)
        end do

    end function set_alpha_beta_spins

    function create_one_spin_basis(n_orbs, n_spins) result(one_spin_basis) 
        integer, intent(in) :: n_orbs, n_spins 
        integer(n_int), allocatable :: one_spin_basis(:) 
#ifdef __DEBUG
        character(*), parameter :: this_routine = "create_one_spin_basis" 
#endif

        integer :: n_max_states, i, right_zero
        integer(n_int) :: count_mask, n_set_zero

        n_max_states = int(choose(n_orbs, n_spins))

        ! 
        allocate(one_spin_basis(n_max_states)) 

        ! create the first basis state: 
        one_spin_basis(1) = maskr(n_spins) 

        ! implement my matlab routine to create the basis states: 
        ! but now i have to do it for spin-orbital encoding!
        ! but do this afterwards so the mapping is easier for off-half-filling!
        do i = 2, n_max_states

            ! copy last state: 
            one_spin_basis(i) = one_spin_basis(i-1)

            ! find the right-most zero with atleast one 1 right of it
            right_zero = right_most_zero(one_spin_basis(i), n_orbs)

            ! if the right-most zero is bigger than n_orbs already, we should 
            ! have been finished already.. 
            ASSERT(right_zero <= n_orbs)

            ! i need to count the number of 1 right of this zero 
            ! so i want a mask where every bit right of right_zero is set 
            ! to one 
            count_mask = maskr(right_zero-1) 
            n_set_zero = popcnt(iand(one_spin_basis(i),count_mask)) 

            ! now i want to set the right_most zero to one  
            one_spin_basis(i) = ibset(one_spin_basis(i), right_zero-1)

            ! and everything to 0 right of it for now 
            one_spin_basis(i) = iand(one_spin_basis(i), not(count_mask)) 

            ! and then we want to set n_set_zero bits to the very right of 
            ! our bit-string 
            one_spin_basis(i) = merge_bits(one_spin_basis(i), maskr(n_set_zero-1,n_int), not(maskr(n_set_zero-1,n_int)))

        end do

    end function create_one_spin_basis

    integer function right_most_zero(i, n_orbs) 
        ! gives the position of the right most zero with atleast one 1 
        ! right of it in an integer bit representation, up to a input 
        ! dependent orbital n_orbs 
        integer(n_int), intent(in) :: i
        integer, intent(in) :: n_orbs
#ifdef __DEBUG
        character(*), parameter :: this_routine = "right_most_zero" 
#endif
        integer, allocatable :: nZeros(:), nOnes(:) 
        integer(n_int) :: j 

        ! first truncate the integer to be save.. 
        ! set all the ones left of n_orbs to 0 
        j = iand(i, maskr(n_orbs))

        ! be sure to avoid the edge case: 
        if (j == 0_n_int) then
            right_most_zero = -1 
            return
        end if

        ! i could misuse the decode_bit_det functionality 
        allocate(nZeros(popcnt(not(j))))
        allocate(nOnes(popcnt(j)))

        ! but exclude the unnecessary 0 left of n_orbs..
        call decode_bit_det(nZeros, [not(j)])
        call decode_bit_det(nOnes, [j])

        ! in this way to encode it.. the right most 0 is then the minumum 
        ! in nZeros which is bigger then the minumum in nOnes
        right_most_zero = minval(nZeros, nZeros > minval(nOnes))

    end function right_most_zero

    function calc_n_double(n_orbs, n_alpha, n_beta) result(n_double)
        ! routine to calculate the number of determinants with different 
        ! number of doubly occupied sites..
        integer, intent(in) :: n_orbs, n_alpha, n_beta
        integer, allocatable :: n_double(:)
#ifdef __DEBUG
        character(*), parameter :: this_routine = "calc_n_double"
#endif

        real(dp) :: n_first
        integer :: n_max, n_min, i
        ! the number of possible open shell determinants is:
        ! for the half-filled ms=0 case, it is how often we can distribute 
        ! n/2 electrons in n orbitals (since the rest then gets filled up 
        ! by the other spin electrons -> nchoosek(n, n/2) 
        ! todo: merge the cc_amplitudes branch or atleast the 
        ! functionality in there! 
        ! for ms/=0 but still at half-filling it i still the binomial 
        ! coefficient of one spin-type 

        ! off-half filling, we get more then one possible distribution of 
        ! the electrons of the other spin for each distribution of the 
        ! first spin type.. 
        ! so we get some product of binomial of the remaining possible 
        ! distributions! 

        ! todo: write those formulas up! 

        ! todo: a way to create it would be like in my matlab code to 
        ! create the heisenberg basis states.. there i even did it in 
        ! a ordered fashion.

        ! the number of ways to distribute alpha spins on the lattice 
!         n_alpha = binomial(nbasis/2, nOccAlpha)

!         num_open_shell_dets = n_alpha * binomial(nbasis/2 - nOccAlpha, nOccBeta) 

        ! maybe the number of determinants with one double occupancy is 
        ! also of interest.. (those are the single excitations from the 
        ! fully open shell dets(in the case off-half-filling it could a
        ! single excitation can also lead to another fully-open-shell det, but
        ! these are already covered.
        ! one beta electron MUST reside in a orbital with a alpha spin (thats 
        ! just nOccAlpha) and the rest (nOccBeta-1) is on an empty orbital
!         num_one_double_occ = n_alpha * nOccAlpha * & 
!             binomial(nbasis/2 - nOccAlpha, nOccBeta - 1)

        ! just for the fun of it: i could calculate the size of each of 
        ! those space(1 double occ, 2 double occ.. etc. 

        ! only do that for less equal half-filling! 
!         ASSERT(n_alpha + n_beta <= n_orbs)

        n_max = max(n_alpha, n_beta) 
        n_min = min(n_alpha, n_beta)

        n_first = choose(n_orbs, n_max)

        allocate(n_double(n_min+1)) 

        do i = 0, n_min 
            n_double(i+1) = int(n_first * choose(n_max, i) * & 
                choose(n_orbs - n_max, n_min - i))
        end do

        ! make sure that the sum of basis states is the whole hilber space
        ASSERT(sum(n_double) == choose(n_orbs, n_alpha)*choose(n_orbs,n_beta))

    end function calc_n_double

    function get_optimal_correlation_factor() result(corr_factor) 
        ! Hongjuns derivation was for the k-space hubbard and in the low 
        ! density and U limit though.. 
        real(dp) :: corr_factor
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "get_optimal_correlation_factor" 
#endif

        ASSERT(associated(lat)) 

        ! the sign is not quite sure here.. which i need to take to 
        ! calculate the hermitian matrix elements..
        corr_factor = -log(abs(real(uhub,dp)/real(4 * lat%get_ndim() * bhub, dp)) + 1.0_dp)


    end function get_optimal_correlation_factor

    subroutine init_tmat(lat)
        class(lattice), optional :: lat

        ! i should create a new, more flexible routine which sets up the 
        ! TMAT for the different lattice types. although i am not sure if 
        ! we need this anymore 
        ! this also depends on the boundary conditions
        character(*), parameter :: this_routine = "init_tmat"

        integer :: i, ind
        ! depending on the input i either create tmat2d here or is have to 
        ! set it up, so it can be used to create the lattice.. 
        ! but for the beginning i think i want to set it up here from the 
        ! lattice structure! 
        ! if the lattice class is alread set up and initialized, this indicates
        ! that it was build-in created and tmat has to be calculated from it 
        ! now!
        if (present(lat)) then 
            ! what do i need to do? 
            ! loop over the indices in the lattice and get the neighbors
            ! and i have to store it in spin-indices remember that!
            if (associated(tmat2d)) deallocate(tmat2d)
            allocate(tmat2d(nbasis,nbasis))
            tmat2d = 0.0_dp

            do i = 1, lat%get_nsites() 
                ind = lat%get_site_index(i)
                associate(next => lat%get_neighbors(i))
                    ! beta orbitals:
                    tmat2d(2*ind - 1, 2*next - 1) = bhub 
                    ! alpha: 
                    tmat2d(2*ind, 2*next) = bhub
                    
                    ASSERT(all(next > 0))
                    ASSERT(all(next <= nbasis/2))
                end associate
                ASSERT(lat%get_nsites() == nbasis/2)
                ASSERT(ind > 0) 
            ASSERT(ind <= nbasis/2)
            
            end do

        else
            ! this indicates that tmat has to be created from an fcidump 
            ! and the lattice is set up afterwards!
        end if

    end subroutine init_tmat
    !
    ! Generic excitaiton generator
    subroutine gen_excit_rs_hubbard (nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                      ex, tParity, pGen, hel, store, run)

        use SystemData, only: nel
        use bit_rep_data, only: NIfTot
        use FciMCData, only: excit_gen_store_type
        use constants, only: n_int, dp, bits_n_int
        use get_excit, only: make_single
        implicit none

        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ic, ex(2,2)
        integer(n_int), intent(out) :: ilutJ(0:NifTot)
        real(dp), intent(out) :: pGen
        logical, intent(out) :: tParity
        HElement_t(dp), intent(out) :: hel
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: run

        character(*), parameter :: this_routine = "gen_excit_rs_hubbard"

        integer :: iunused, ind , elec, id, src, orb, n_avail, n_orbs, i
        integer, allocatable :: neighbors(:), orbs(:)
        real(dp), allocatable :: cum_arr(:)
        real(dp) :: cum_sum, elem, r, p_elec, p_orb

        iunused = exflag; 

        ASSERT(associated(lat))

        ic = 1
        ! i only have single excitations in the hubbard model 
        ! the first plan is to choose an electron at random 
        elec = 1 + int(genrand_real2_dsfmt() * nel) 

        p_elec = 1.0_dp / real(nel, dp)
        ! and then from the neighbors of this electron we pick an empty 
        ! spinorbital randomly, since all have the same matrix element 
        src = nI(elec) 

        ! get the spatial index 
        id = get_spatial(src)
!         id = gtid(src) 

        ! now get neighbors
!         n_orbs = lat%get_num_neighbors(id)
        neighbors = lat%get_spinorb_neighbors(src)

!         allocate(neighbors(n_orbs)) 
!         allocate(orbs(n_orbs))
!         neighbors = lat%get_neighbors(id) 

!         if (is_beta(src)) then 
!             neighbors = 2 * neighbors - 1
!         else 
!             neighbors = 2 * neighbors
!         end if

!         if (t_trans_corr) then 
            call create_cum_list_rs_hubbard(ilutI, src, neighbors, cum_arr, cum_sum)

            if (cum_sum < EPS) then 
                nJ(1) = 0
                pgen = 0.0_dp 
                return
            end if
            r = genrand_real2_dsfmt() * cum_sum 
            ind = binary_search_first_ge(cum_arr, r)

            if (ind == 1) then 
                p_orb = cum_arr(1) / cum_sum
            else
                p_orb = (cum_arr(ind) - cum_arr(ind-1)) / cum_sum 
            end if

            orb = neighbors(ind)
!         else
! 
!             call create_avail_neighbors_list(ilutI, neighbors, orbs, n_avail)
! 
!             if (n_avail == 0) then 
!                 nJ(1) = 0
!                 pgen = 0.0_dp 
!                 return
!             end if
! 
!             ind = 1 + int(genrand_real2_dsfmt() * n_avail)
!             p_orb = 1.0_dp / real(n_avail, dp) 
! 
!             orb = orbs(ind) 
! 
!         end if

        pgen = p_elec * p_orb

        call make_single(nI, nJ, elec, orb, ex, tParity) 

        ilutJ = ilutI 
        clr_orb(ilutJ, src)
        set_orb(ilutJ, orb)

    end subroutine gen_excit_rs_hubbard

    function calc_pgen_rs_hubbard(nI, ilutI, ex, ic) result(pgen) 
        ! i also need a pgen recalculator.. specifically for the HPHF 
        ! implementation and i need to take the transcorrelated keyword 
        ! into account here! 
        integer, intent(in) :: nI(nel), ex(2,2), ic 
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        real(dp) :: pgen 
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "calc_pgen_rs_hubbard"
#endif
        integer :: src, tgt, n_orbs
        real(dp) :: p_elec, p_orb, cum_sum
        real(dp), allocatable :: cum_arr(:)
        integer, allocatable :: orbs(:)

        ! only single excitations in the real-space hubbard
        if (ic /= 1) then 
            pgen = 0.0_dp 
            return
        end if

        src = ex(1,1)
        tgt = ex(2,1)

        ! can i assert the same spin of the 2 involved orbitals? 
        ! just return 0 if both have different spin? 
        ASSERT(is_beta(src) .eqv. is_beta(tgt))
        ! and assert that we actually take a valid excitation:
        ASSERT(any(tgt == lat%get_spinorb_neighbors(src)))
        ASSERT(IsOcc(ilutI, src))
        ASSERT(IsNotOcc(ilutI, tgt))

        p_elec = 1.0_dp / real(nel, dp) 

!         if (t_trans_corr) then 
            call create_cum_list_rs_hubbard(ilutI, src, lat%get_spinorb_neighbors(src), &
                cum_arr, cum_sum, tgt, p_orb) 
            if (cum_sum < EPS) then 
                pgen = 0.0_dp 
                return
            end if

!         else 
!             ! i should also write a routine which gives me the 
!             ! neighboring orbitals and the number of possible hops 
!             call create_avail_neighbors_list(ilutI, lat%get_spinorb_neighbors(src), & 
!                 orbs, n_orbs) 
! 
!                 p_orb = 1.0_dp / real(n_orbs, dp) 
! 
!         end if

        pgen = p_elec * p_orb

    end function calc_pgen_rs_hubbard

    subroutine create_cum_list_rs_hubbard(ilutI, src, neighbors, cum_arr, cum_sum, &
            tgt, cpt)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(in) :: src, neighbors(:)
        real(dp), intent(out) :: cum_sum 
        real(dp), intent(out), allocatable :: cum_arr(:)
        integer, intent(in), optional :: tgt
        real(dp), intent(out), optional :: cpt
#ifdef __DEBUG
        character(*), parameter :: this_routine = "create_cum_list_rs_hubbard"
#endif

        real(dp) :: elem
        integer :: i, nI(nel)

        ASSERT(IsOcc(ilutI,src))

        call decode_bit_det(nI, ilutI)

        allocate(cum_arr(size(neighbors)))
        cum_arr = 0.0_dp
        cum_sum = 0.0_dp
        if (present(tgt)) then 
            ASSERT(present(cpt))
            do i = 1, ubound(neighbors,1)
                elem = 0.0_dp
                ASSERT(is_beta(src) .eqv. is_beta(neighbors(i)))
                if (IsNotOcc(ilutI, neighbors(i))) then 
                    elem = abs(get_offdiag_helement_rs_hub(nI,[src,neighbors(i)],.false.))
                end if
                cum_sum = cum_sum + elem
                if (neighbors(i) == tgt) then 
                    cpt = elem
                end if
            end do
            if (cum_sum < EPS) then
                cpt = 0.0
            else 
                cpt = cpt / cum_sum
            end if
        else
            do i = 1, ubound(neighbors,1)
                elem = 0.0_dp
                ASSERT(is_beta(src) .eqv. is_beta(neighbors(i)))
                if (IsNotOcc(ilutI,neighbors(i))) then 
                    elem = abs(get_offdiag_helement_rs_hub(nI,[src,neighbors(i)],.false.))
                end if
                cum_sum = cum_sum + elem
                cum_arr(i) = cum_sum
            end do
        end if
        
    end subroutine create_cum_list_rs_hubbard

    subroutine create_avail_neighbors_list(ilutI, neighbors, orbs, n_orbs)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(in) :: neighbors(:) 
        integer, intent(out), allocatable :: orbs(:)
        integer, intent(out) :: n_orbs 

        integer :: i, temp_orbs(size(neighbors))

        n_orbs = 0
        temp_orbs = 0

        do i = 1, ubound(neighbors,1)
            if (IsNotOcc(ilutI,neighbors(i))) then 
                n_orbs = n_orbs + 1
                temp_orbs(n_orbs) = neighbors(i)
            end if
        end do

        allocate(orbs(n_orbs), source=temp_orbs(1:n_orbs))

    end subroutine create_avail_neighbors_list

    function trans_corr_fac(ilutI, src, tgt) result(weight)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(in) :: src, tgt
        real(dp) :: weight 
#ifdef __DEBUG
        character(*), parameter :: this_routine = "trans_corr_fac"
#endif
        real(dp) :: ni_opp, nj_opp 

        ! if the spins are not the same, something went wrong..
        ASSERT(is_beta(src) .eqv. is_beta(tgt))

        ni_opp = 0.0_dp
        nj_opp = 0.0_dp

        if (is_beta(src)) then 
            ! check if alpha orbital (i) and (j) in ilutI is occupied
            if (IsOcc(ilutI,get_alpha(src))) then 
                nj_opp = 1.0_dp
            end if
            if (IsOcc(ilutI,get_alpha(tgt))) ni_opp = 1.0_dp
        else 
            if (IsOcc(ilutI,get_beta(src))) nj_opp = 1.0_dp
            if (IsOcc(ilutI,get_beta(tgt))) ni_opp = 1.0_dp
        end if

        weight = exp(trans_corr_param*(nj_opp - ni_opp))

    end function trans_corr_fac

    function get_helement_rs_hub_ex_mat(nI, ic, ex, tpar) result(hel)
        integer, intent(in) :: nI(nel)
        integer, intent(in) :: ic, ex(2,ic)
        logical, intent(in) :: tpar
        HElement_t(dp) :: hel 

        if (ic == 0) then 
            ! diagonal matrix element -> sum over doubly occupied orbitals! 
            hel = get_diag_helemen_rs_hub(nI) 

        else if (ic == 1) then 
            ! one-body operator: 
            ! here we need to make the distinction, if we are doing a 
            ! transcorrelated hamiltonian or not 
            hel = get_offdiag_helement_rs_hub(nI, ex(:,1), tpar)

        else 
            ! zero matrix element! 
            hel = h_cast(0.0_dp)

        end if

    end function get_helement_rs_hub_ex_mat

    function get_helement_rs_hub_general(nI, nJ, ic_ret) result(hel)
        integer, intent(in) :: nI(nel), nJ(nel) 
        integer, intent(inout), optional :: ic_ret 
        HElement_t(dp) :: hel
        
        integer :: ic, ex(2,2)
        logical :: tpar
        integer(n_int) :: ilutI(0:NIfTot), ilutJ(0:NIfTot)

        if (present(ic_ret)) then 
            if (ic_ret == 0) then 
                hel = get_diag_helemen_rs_hub(nI)

            else if (ic_ret == 1) then 
                ex(1,1) = 1
                call GetExcitation(nI, nJ, nel, ex, tpar)
                hel = get_offdiag_helement_rs_hub(nI, ex(:,1), tpar) 

            else if (ic_ret == -1) then 
                ! this indicates that ic_ret wants to get returned instead of 
                ! beeing calculated 
                ! its the same as if no ic_ret is present 
                call EncodeBitDet(nI, ilutI)
                call EncodeBitDet(nJ, ilutJ)

                ic_ret = FindBitExcitLevel(ilutI, ilutJ)

                if (ic_ret == 0) then 
                    hel = get_diag_helemen_rs_hub(nI)
                else if (ic_ret == 1) then 
                    ex(1,1) = 1
                    call GetBitExcitation(ilutI, ilutJ, ex, tpar)

                    hel = get_offdiag_helement_rs_hub(nI, ex(:,1), tpar)

                else 
                    hel = h_cast(0.0_dp)

                end if
            else
                hel = h_cast(0.0_dp)
            end if
        else
            call EncodeBitDet(nI, ilutI)
            call EncodeBitDet(nJ, ilutJ)

            ic = FindBitExcitLevel(ilutI, ilutJ)

            if (ic == 0) then 
                hel = get_diag_helemen_rs_hub(nI)
            else if (ic == 1) then 
                ex(1,1) = 1
                call GetBitExcitation(ilutI, ilutJ, ex, tpar)
                hel = get_offdiag_helement_rs_hub(nI, ex(:,1), tpar) 

            else
                hel = h_cast(0.0_dp)
            end if
        end if
    
    end function get_helement_rs_hub_general

    ! also optimize the matrix element calculation
    function get_diag_helemen_rs_hub(nI) result(hel)
        use double_occ_mod, only: count_double_orbs
        integer, intent(in) :: nI(nel)
        HElement_t(dp) :: hel

        integer(n_int) :: ilut(0:NIfTot)
        ! the diagonal matrix element is essentialy just the number of 
        ! doubly occupied sites times U

        call EncodeBitDet(nI, ilut)

        hel = h_cast(uhub * count_double_orbs(ilut))

    end function get_diag_helemen_rs_hub

    function get_offdiag_helement_rs_hub(nI, ex, tpar) result(hel)
        integer, intent(in) :: nI(nel), ex(2)
        logical, intent(in) :: tpar
        HElement_t(dp) :: hel

        integer(n_int) :: ilut(0:NIfTot)
        real(dp) :: n_i, n_j

        ! in case we need it, the off-diagonal, except parity is just 
        ! -t if the hop is possible
        hel = GetTMatEl(ex(1),ex(2))

        if (tpar) hel = -hel

        ! put the transcorrelated stuff here for now, alhtough it would be 
        ! better to do it as a procedure pointer.. 
        if (t_trans_corr) then 
            call EncodeBitDet(nI, ilut)
            if (t_trans_corr_new) then
                n_j = get_opp_spin(ilut, ex(1))
                n_i = get_opp_spin(ilut, ex(2))

                ! try to go one step further before the transformation to 
                ! the momentum space.. and check if the results are still 
                ! correct.. 

                hel = hel * (1.0_dp + (exp(trans_corr_param) - 1.0_dp) * n_j + & 
                                      (exp(-trans_corr_param)- 1.0_dp) * n_i - & 
                                       2.0_dp*(cosh(trans_corr_param) - 1.0_dp)*n_i*n_j)

           else

                hel = hel * exp(trans_corr_param * & 
                    (get_opp_spin(ilut, ex(1)) - get_opp_spin(ilut, ex(2))))
            end if

        end if

        ! it is not really a 2-body trans, but still use the flag and 
        ! parameter
        if (t_trans_corr_2body) then 
            call EncodeBitDet(nI, ilut) 
            hel = hel * exp(trans_corr_param_2body * & 
                (get_spin_opp_neighbors(ilut, ex(1)) - get_spin_opp_neighbors(ilut,ex(2))))
        end if

    end function get_offdiag_helement_rs_hub

    function get_opp_spin(ilut, spin_orb) result(opp_spin) 
        integer(n_int), intent(in) :: ilut(0:NIfTot) 
        integer, intent(in) :: spin_orb
        real(dp) :: opp_spin 

        opp_spin = 0.0_dp

        if (is_beta(spin_orb)) then 
            if (IsOcc(ilut, spin_orb + 1)) opp_spin = 1.0_dp 
        else 
            if (IsOcc(ilut, spin_orb -1 )) opp_spin = 1.0_dp
        end if

    end function get_opp_spin

    function get_spin_opp_neighbors(ilut, spin_orb) result(spin_opp_neighbors) 
        ! function to give the number of opposite spin electron neighbors 
        integer(n_int), intent(in) :: ilut(0:niftot) 
        integer, intent(in) :: spin_orb 
        real(dp) :: spin_opp_neighbors
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "get_spin_opp_neighbors" 
#endif 
        integer :: i, flip 
        integer, allocatable :: neighbors(:) 

        ASSERT(associated(lat))

        spin_opp_neighbors = 0.0_dp 


        ! get the spin-opposite neigbhors 
        if (is_beta(spin_orb)) then 
            neighbors = lat%get_spinorb_neighbors(spin_orb) + 1
        else 
            neighbors = lat%get_spinorb_neighbors(spin_orb) - 1
        end if

        do i = 1, size(neighbors) 
            if (IsOcc(ilut, neighbors(i))) spin_opp_neighbors = spin_opp_neighbors + 1.0_dp
        end do

    end function get_spin_opp_neighbors


    ! what else?
    function create_neel_state(ilut_neel) result(neel_state)
        ! probably a good idea to have a routine which creates a neel state
        ! (or one of them if it is ambigous)
        integer(n_int), intent(out), optional :: ilut_neel(0:NIfTot)
        integer :: neel_state(nel)
        character(*), parameter :: this_routine = "create_neel_state"

        integer :: i, j, k, l, spin, ind
        

        if (tHPHF) then 
            call stop_all(this_routine, "remember to change for HPHF!")
        end if
        ! make this independent of the lattice mod, to call it as early 
        ! as possible in the initialization
!         ASSERT(associated(lat))

        ! also assert that we have at most half-filling.. otherwise it does 
        ! not make so much sense to do this 
!         ASSERT(nel <= nBasis/2)
!         ASSERT(nel <= lat%get_nsites())
        ! there is no neel state for all lattice, but atleast something 
        ! similar to it.. 
        ! atleast for the lattice where there is a neel state, i should 
        ! create it automaticall and for the other a similar one 
        ! the lattice should already be intialized 
        ! i was thinking of putting this functionality into the 
        ! lattice_mod, but i guess i should keep it seperate, since it 
        ! actually has to do more with the system or?
        neel_state = 0

        select case (lattice_type) 

        case ('chain')
            neel_state = create_neel_state_chain()

        case ('square','rectangle','triangle','triangular','hexagonal','kagome')
            ! check if length_x is mod 2 
            if (mod(length_x, 2) == 0) then 

                neel_state = create_neel_state_chain()

                ! and flip the spins in every other column 
                do i = 1, length_y, 2
                    ! loop over every second column 
                    if (i*length_x >= nel) exit
                    do j = 1, length_x
                        ! whats the total index? 
                        ind = i*length_x + j
                        if (ind > nel) exit

                        ! flip the spins.. this should be way easier.. 
                        if (is_beta(neel_state(ind))) then
                            neel_state(ind) = neel_state(ind) + 1
                        else
                            neel_state(ind) = neel_state(ind) - 1
                        end if
                    end do
                end do

            else 
                ! here it is easy it is just like the chain case 
                neel_state = create_neel_state_chain()

            end if

        case ('cube') 
            ! not yet implemented
            ASSERT(.false.)

        case ('tilted','tilted-square','square-tilted')
            ! do is similar to the actual construction of the lattice 
            ! but way nicer as this stuff below..
            k = 0
            l = 1 
            spin = 1
            do i = -length_x+1, 0
                do j = -k, k

                    neel_state(l) = spin 
                    l = l + 1
                    if (l > nel) exit 

                    if (is_beta(spin)) then 
                        spin = spin + 3 
                    else 
                        spin = spin + 1 
                    end if 
                end do
                k = k + 1 
                if (l > nel) exit 

                ! in the first half we need the same starting spin 
                spin = spin - 1

            end do
            k = k - 1

            ! finished already?
            if (l > nel) return

            spin = spin + 1

            do i = 1, length_x 
                do j = -k, k

                    neel_state(l) = spin 

                    l = l + 1 

                    if (l > nel) exit 

                    if (is_beta(spin)) then 
                        spin = spin + 3 
                    else 
                        spin = spin + 1
                    end if
                end do
                k = k - 1 
                if (l > nel) exit 

                spin = spin + 1

            end do

        case default 
            call stop_all(this_routine, "unknown lattice type!")
        end select 

        if (present(ilut_neel)) then 
            call EncodeBitDet(neel_state, ilut_neel)
        end if
    end function create_neel_state

    function create_neel_state_chain() result(neel_state)
        integer :: neel_state(nel)

        integer :: i

        neel_state = [(i, i = 1, 2*nel-1,2)]
        neel_state(2:nel:2) = neel_state(2:nel:2) + 1

    end function create_neel_state_chain

    function get_umat_el_hub(i,j,k,l) result(hel)
        integer, intent(in) :: i, j, k, l
        HElement_t(dp) :: hel 
#ifdef __DEBUG
        character(*), parameter :: this_routine = "get_umat_el_hub"
#endif
        if (i == j .and. i == k .and. i == l) then 
            hel = h_cast(uhub)
        else 
            hel = h_cast(0.0_dp)
        end if

        ASSERT(i > 0)
        ASSERT(i <= nbasis/2)
        ASSERT(j > 0) 
        ASSERT(j <= nbasis/2)
        ASSERT(k > 0) 
        ASSERT(k <= nbasis/2)
        ASSERT(l > 0) 
        ASSERT(l <= nbasis/2)

    end function get_umat_el_hub

end module real_space_hubbard
        
        
