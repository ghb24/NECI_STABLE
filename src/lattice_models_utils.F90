#include "macros.h"

! create a module for useful small functions used in various lattice 
! model implementations and modules to avoid circular dependencies 

module lattice_models_utils 

    use constants, only: dp, n_int, bits_n_int, eps, pi, lenof_sign

    use util_mod, only: binary_search, binary_search_first_ge, choose, swap

    use sort_mod, only: sort

    use lattice_mod, only: lat, sort_unique, get_helement_lattice

    use SystemData, only: symmetry, G1, nel, nbasis, nBasisMax, t_k_space_hubbard, &
                          tHub, nOccBeta, nOccAlpha, t_trans_corr_2body, tHPHF, &
                          lattice_type, length_x, length_y, t_trans_corr_hop, &
                          arr, brr

    use symdata, only: symtable, SymConjTab

    use bit_rep_data, only: niftot, nifd

    use DetBitOps, only: ilut_lt, ilut_gt, EncodeBitDet, TestClosedShellDet, & 
                         DetBitLT, MaskAlpha, MaskBeta

    use bit_reps, only: decode_bit_det, extract_sign

    use SymExcitDataMod, only: KPointToBasisFn

    use sym_mod, only: mompbcsym

    use dSFMT_interface, only: genrand_real2_dsfmt

    use Parallel_neci, only: iProcIndex, root

    use get_excit, only: make_double

    implicit none

    interface swap_excitations 
        module procedure swap_excitations_higher
        module procedure swap_excitations_singles
    end interface swap_excitations

contains 


    function return_hphf_sym_det(ilut_in) result(ilut_out)
        ! to avoid circular dependencies and due to the strange implementation
        ! to find the symmetry conjugated determinant of an HPHF pair 
        ! create a new routine to return a open-shell determinant where the 
        ! last single occupied spatial orbital is an alpha spin
        ! this is the convention in the storage of the hphfs
        ! this can easily be tested by checking if the bit-encoded determinant 
        ! has an higher integer value! 
        ! different to the original implementation of this routine 
        ! standardyl we only return the determinant which should be stored 
        ! in the CurrentDets. so if ilut_in is already this determinant 
        ! ilut_out will be == ilut_in
        ! and it also deals with closed-shell dets, where it will just 
        ! return the same determinant 
        integer(n_int), intent(in) :: ilut_in(0:niftot)
        integer(n_int) :: ilut_out(0:niftot)
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "return_hphf_sym_det"
#endif
        INTEGER(n_int) :: iLutAlpha(0:NIfTot),iLutBeta(0:NIfTot)
        INTEGER :: i

        if (TestClosedShellDet(ilut_in)) then 
            ilut_out = ilut_in 
            return 
        end if

        ilut_out(:)=0_n_int
        iLutAlpha(:)=0_n_int
        iLutBeta(:)=0_n_int

        ! this is taken from HPHFRandExcitMod
        do i=0,NIfD
            !Seperate the alpha and beta bit strings
            iLutAlpha(i)=IAND(ilut_in(i),MaskAlpha)    
            iLutBeta(i)=IAND(ilut_in(i),MaskBeta)

            !Shift all alpha bits to the left by one.
            iLutAlpha(i)=ISHFT(iLutAlpha(i),-1)  
            !Shift all beta bits to the right by one.
            iLutBeta(i)=ISHFT(iLutBeta(i),1)   
            !Combine the bit strings to give the final bit representation.
            ilut_out(i)=IOR(iLutAlpha(i),iLutBeta(i))    
        end do

        i = DetBitLT(ilut_in, ilut_out)

        ! i == 1 indicated that ilut_in is "less" than the symmetric
        ! so ilut_out is the to be stored one 
        if (i == -1) then 
            ilut_out = ilut_in
        end if

    end function return_hphf_sym_det

    subroutine swap_excitations_higher(nI, ex, nJ, ex2) 
        ! routine to quickly, without make_double and make_triple 
        ! create the excited determinant nJ and ex2 to go from nJ -> nI 
        ! for the test of the order of matrix element calculation in the 
        ! transcorrelated approach. due to non-hermiticity 
        ! <i|H|j> /= <j|H|i> 
        ! also make this work for triple excitations! 
        integer, intent(in) :: nI(nel), ex(:,:) 
        integer, intent(out) :: nJ(nel) 
        integer, intent(out), allocatable :: ex2(:,:)
#ifdef __DEBUG
        character(*), parameter :: this_routine = "swap_excitations" 
#endif

        ASSERT(size(ex,1) == 2) 
        ASSERT(size(ex,2) == 2 .or. size(ex,2) == 3)

        allocate(ex2(size(ex,1),size(ex,2)), source = ex)

        nJ = nI
        where (nJ == ex(1,1)) nJ = ex(2,1)
        where (nJ == ex(1,2)) nJ = ex(2,2)
        if (size(ex,2) == 3) then 
            where (nJ == ex(1,3)) nJ = ex(2,3)
        end if

        call sort(nJ)
        ex2 = ex 
        call swap(ex2(1,1),ex2(2,1))
        call swap(ex2(1,2),ex2(2,2))
        if (size(ex,2) == 3) then 
            call swap(ex2(1,3),ex2(2,3))
        end if

    end subroutine swap_excitations_higher

    subroutine swap_excitations_singles(nI, ex, nJ, ex2)
        integer, intent(in) :: nI(nel), ex(2) 
        integer, intent(out) :: nJ(nel), ex2(2)
        
        nJ = nI 

        where (nJ == ex(1)) nJ = ex(2)

        call sort(nJ) 

        ex2 = ex

        call swap(ex2(1),ex2(2))

    end subroutine swap_excitations_singles

    subroutine pick_spin_opp_elecs(nI, elecs, p_elec) 
        integer, intent(in) :: nI(nel)
        integer, intent(out) :: elecs(2)
        real(dp), intent(out) :: p_elec
#ifdef __DEBUG
        character(*), parameter :: this_routine = "pick_spin_opp_elecs"
#endif
        ! think of a routine to get the possible spin-opposite electron 
        ! pairs. i think i could do that way more efficiently, but do it in 
        ! the simple loop way for now 
        do 
            elecs(1) = 1 + int(genrand_real2_dsfmt() * nel)
            do 
                elecs(2) = 1 + int(genrand_real2_dsfmt() * nel) 

                if (elecs(1) /= elecs(2)) exit 

            end do
            if (get_ispn(nI(elecs)) == 2) exit
        end do

        ! output in ordered form 
        elecs = [minval(elecs),maxval(elecs)]

        ! actually the probability is twice that or? 
        ! or doesnt that matter, since it is the same
        p_elec = 1.0_dp / real(nOccBeta * nOccAlpha, dp)

    end subroutine pick_spin_opp_elecs

    subroutine pick_from_cum_list(cum_arr, cum_sum, ind, pgen) 
        real(dp), intent(in) :: cum_arr(:), cum_sum
        integer, intent(out) :: ind
        real(dp), intent(out) :: pgen 
#ifdef __DEBUG
        character(*), parameter :: this_routine = "pick_from_cum_list" 
#endif
        real(dp) :: r

        if (cum_sum < EPS) then 
            ind = -1 
            pgen = 0.0_dp
            return 
        end if

        r = genrand_real2_dsfmt() * cum_sum

        ind = binary_search_first_ge(cum_arr, r) 

        if (ind == 1) then 
            pgen = cum_arr(1)/cum_sum 
        else 
            pgen = (cum_arr(ind) - cum_arr(ind - 1)) / cum_sum
        end if

    end subroutine pick_from_cum_list

    subroutine gen_all_excits_r_space_hubbard(nI, n_excits, det_list) 
        ! for the purpose of excitation generation and time-step and 
        ! pDoubles determination create a routine to create all possible 
        ! excitations for a given determinant nI 
        ! the system specific variables need to be determined already 
        ! before! 
        integer, intent(in) :: nI(nel)
        integer, intent(out) :: n_excits
        integer(n_int), intent(out), allocatable :: det_list(:,:) 
        character(*), parameter :: this_routine = "gen_all_excits_r_space_hubbard"

        integer :: save_excits, n_doubles, i
        integer(n_int), allocatable :: double_dets(:,:), temp_dets(:,:)

        call gen_all_singles_rs_hub(nI, n_excits, det_list) 

        ! if we have hopping transcorr we also have to compute the 
        ! possible doubles 
        ! we know does are different than any single
        if (t_trans_corr_hop) then 
            save_excits = n_excits

            call gen_all_doubles_rs_hub_hop_transcorr(nI, n_doubles, double_dets) 

            n_excits = n_excits + n_doubles

            allocate(temp_dets(0:niftot, save_excits), source = det_list(:,1:save_excits))

            deallocate(det_list) 

            allocate(det_list(0:niftot,n_excits)) 

            det_list(:,1:save_excits) = temp_dets 

            det_list(:,save_excits+1:n_excits) = double_dets

        end if

        call sort(det_list, ilut_lt, ilut_gt)

        ! if we have HPHF turned on we want to "spin-purify" the excitation 
        ! list 
!         print *, "un-purified basis: "
!         do i = 1, n_excits
!             call writebitdet(6, det_list(:,i), .true.)
!         end do
        if (tHPHF) then 
            save_excits = n_excits

            if (allocated(temp_dets)) deallocate(temp_dets)

            allocate(temp_dets(0:NIfTot, size(det_list,2)), source = det_list)

            call spin_purify(save_excits, temp_dets, n_excits, det_list)

        end if

!         print *, "purified basis: "
!         do i = 1, n_excits
!             call writebitdet(6, det_list(:,i), .true.)
!         end do
    end subroutine gen_all_excits_r_space_hubbard

    subroutine spin_purify(n_excits_in, det_list_in, n_excits_out, det_list_out)
        ! routine to remove determinants, belonging to the same 
        ! coupled HPHF function 
        integer, intent(in) :: n_excits_in
        integer(n_int), intent(in) :: det_list_in(0:NIfTot,n_excits_in)
        integer, intent(out) :: n_excits_out
        integer(n_int), intent(out), allocatable :: det_list_out(:,:)
#ifdef __DEBUG
        character(*), parameter :: this_routine = "spin_purify"
#endif 
        integer :: nI(nel), nJ(nel), i, pos, cnt
        integer(n_int) :: ilut(0:NIfTot), ilut_sym(0:NIfTot)
        logical :: t_swapped
        integer(n_int), allocatable :: temp_dets(:,:)

        allocate(temp_dets(0:NIfTot,n_excits_in))
        temp_dets = 0_n_int

        cnt = 0

        do i = 1, n_excits_in
        
            ilut = det_list_in(:,i) 

            ! ilut will always be the one, which should be stored in the 
            ! det-list, if i undertand the function below correctly
            ilut_sym = return_hphf_sym_det(ilut)

            pos = binary_search(temp_dets(:,1:cnt), ilut, nifd + 1)

            if (pos < 0) then 
                ! then we have to store it 
                cnt = cnt + 1
                temp_dets(:,cnt) = ilut_sym
                call sort(temp_dets(:,1:cnt), ilut_lt, ilut_gt)
            end if
        end do

        n_excits_out = cnt
        allocate(det_list_out(0:NIfTot,n_excits_out), source = temp_dets(:,1:n_excits_out))

        call sort(det_list_out, ilut_lt, ilut_gt)

    end subroutine spin_purify

    subroutine gen_all_doubles_rs_hub_hop_transcorr(nI, n_excits, det_list) 
        integer, intent(in) :: nI(nel)
        integer, intent(out) :: n_excits
        integer(n_int), intent(out), allocatable :: det_list(:,:) 
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "gen_all_doubles_rs_hub_hop_transcorr"
#endif
        integer(n_int) :: ilut(0:NIfTot), ilutJ(0:NIfTot)
        integer :: n_bound, i, src(2), j, neigh, ex(2,2), pos, a, b
        integer(n_int), allocatable :: temp_list(:,:)
        integer, allocatable :: neighbors(:)
        real(dp) :: elem

        call EncodeBitDet(nI, ilut) 

        n_excits = 1
        
        n_bound = nel*(nel-1)*(nbasis-nel)*(nbasis-nel-1)

        allocate(temp_list(0:NIfTot,n_bound))
        temp_list = 0_n_int

        ! we just have opposite spin doubles! 
        do i = 1, nel - 1
            do j = i + 1, nel 

                src = nI([i,j])

                if (same_spin(src(1),src(2))) cycle 

                ex(1,:) = src 
                ! and now loop over all possible empty spin opposite 
                ! orbitals in a unique way.. 
                do a = 1, nbasis - 1
                    if (IsOcc(ilut,a)) cycle 
                    do b = a + 1, nbasis 
                        if (IsOcc(ilut,b) .or. same_spin(a,b)) cycle 

                        ex(2,:) = [a,b] 

                        elem = abs(get_helement_lattice(nI, 2, ex, .false.))
                        if (elem > EPS) then
!                         if (abs(get_double_helem_rs_hub_transcorr(nI, ex, .false.)) > EPS) then

                            ilutJ = make_ilutJ(ilut, ex, 2)
                            ! actually a search is not really necessary.. since 
                            ! all the single excitations are unique.. but 
                            ! just to be sure
                            pos = binary_search(temp_list(:,1:n_excits), ilutJ, nifd+1)

                            if (pos < 0) then 

                                temp_list(:,n_excits) = ilutJ
                                call sort(temp_list(:,1:n_excits),ilut_lt,ilut_gt)
                                n_excits = n_excits + 1
                                ! damn.. i have to sort everytime i guess..
                            end if
                        end if
                    end do
                end do
            end do
        end do

        n_excits = n_excits - 1
        allocate(det_list(0:NIfTot,n_excits), source = temp_list(:,1:n_excits))

        call sort(det_list, ilut_lt, ilut_gt)

    end subroutine gen_all_doubles_rs_hub_hop_transcorr

    subroutine gen_all_singles_rs_hub(nI, n_excits, det_list)
        ! create all single excitations in the real-space hubbard 
        ! without hopping transcorrelation this is quite easy.. 
        integer, intent(in) :: nI(nel)
        integer, intent(out) :: n_excits
        integer(n_int), intent(out), allocatable :: det_list(:,:) 
        character(*), parameter :: this_routine = "gen_all_singles_rs_hub" 

        if (t_trans_corr_hop) then 
            call gen_all_singles_rs_hub_hop_transcorr(nI, n_excits, det_list)
        else 
            call gen_all_singles_rs_hub_default(nI, n_excits, det_list) 
        end if

    end subroutine gen_all_singles_rs_hub

    subroutine gen_all_singles_rs_hub_hop_transcorr(nI, n_excits, det_list)
        integer, intent(in) :: nI(nel)
        integer, intent(out) :: n_excits
        integer(n_int), intent(out), allocatable :: det_list(:,:) 
#ifdef __DEBUG
        character(*), parameter :: this_routine = "gen_all_singles_rs_hub_hop_transcorr"
#endif 
        integer(n_int) :: ilut(0:NIfTot), ilutJ(0:NIfTot)
        integer :: n_bound, i, src, j, neigh, ex(2,2), pos, a
        integer(n_int), allocatable :: temp_list(:,:)
        integer, allocatable :: neighbors(:)
        real(dp) :: elem

        call EncodeBitDet(nI, ilut) 

        n_excits = 1

        n_bound = nel * (nbasis - nel)
        allocate(temp_list(0:NIfTot,n_bound))
        temp_list = 0_n_int
        ex = 0

        ! loop over electrons 
        do i = 1, nel 
            ! but now loop over all orbitals 
            src = nI(i)
            ex(1,1) = src 

            do a = 1, nbasis

                ! only same-spin excitations possible
                if (same_spin(src,a) .and. IsNotOcc(ilut,a)) then 
                    ex(2,1) = a
                    elem = abs(get_helement_lattice(nI, 1, ex, .false.))
!                     elem = abs(get_single_helem_rs_hub_transcorr(nI, ex, .false.))

                    if (elem > EPS) then 

                        ilutJ = make_ilutJ(ilut, ex, 1)

                        ! actually a search is not really necessary.. since 
                        ! all the single excitations are unique.. but 
                        ! just to be sure
                        pos = binary_search(temp_list(:,1:n_excits), ilutJ, nifd+1)

                        if (pos < 0) then 

                            temp_list(:,n_excits) = ilutJ
                            call sort(temp_list(:,1:n_excits),ilut_lt,ilut_gt)
                            n_excits = n_excits + 1
                            ! damn.. i have to sort everytime i guess..
                        end if
                    end if
                end if
            end do
        end do

        n_excits = n_excits - 1
        allocate(det_list(0:NIfTot,n_excits), source = temp_list(:,1:n_excits))

        call sort(det_list, ilut_lt, ilut_gt)

    end subroutine gen_all_singles_rs_hub_hop_transcorr

    subroutine gen_all_singles_rs_hub_default(nI, n_excits, det_list) 
        integer, intent(in) :: nI(nel)
        integer, intent(out) :: n_excits
        integer(n_int), intent(out), allocatable :: det_list(:,:) 
#ifdef __DEBUG
        character(*), parameter :: this_routine = "gen_all_singles_rs_hub_default" 
#endif
        integer(n_int) :: ilut(0:NIfTot), ilutJ(0:NIfTot)
        integer :: n_bound, i, src, j, neigh, ex(2,2), pos
        integer(n_int), allocatable :: temp_list(:,:)
        integer, allocatable :: neighbors(:)
        real(dp) :: elem

        ASSERT(associated(lat))

        call EncodeBitDet(nI, ilut) 

        n_excits = 1

        n_bound = nel * (nbasis - nel)
        allocate(temp_list(0:NIfTot,n_bound))
        temp_list = 0_n_int

        ex = 0
        ! loop over electrons
        do i = 1, nel 
            src = nI(i)

            ex(1,1) = src
            neighbors = lat%get_spinorb_neighbors(src)

            ! and loop over the neighboring sites of this electron
            do j = 1, size(neighbors)

                neigh = neighbors(j)

                ex(2,1) = neigh

                ASSERT(same_spin(src,neigh))

                ! if it is not occupied it should be a possible excitation
                if (IsNotOcc(ilut,neighbors(j))) then 
                    ! but to be sure, check the matrix element: 
                    ! but use the lattice get_helement_lattice to 
                    ! avoid circular dependencies
!                     elem = abs(get_offdiag_helement_rs_hub(nI,ex,.false.))
                    elem = abs(get_helement_lattice(nI, 1, ex, .false.))

                    if (elem > EPS) then 

                        ilutJ = make_ilutJ(ilut, ex, 1)

                        ! actually a search is not really necessary.. since 
                        ! all the single excitations are unique.. but 
                        ! just to be sure
                        pos = binary_search(temp_list(:,1:n_excits), ilutJ, nifd+1)

                        if (pos < 0) then 

                            temp_list(:,n_excits) = ilutJ
                            call sort(temp_list(:,1:n_excits),ilut_lt,ilut_gt)
                            n_excits = n_excits + 1
                            ! damn.. i have to sort everytime i guess..
                        end if
                    end if
                end if
            end do
        end do

        n_excits = n_excits - 1
        allocate(det_list(0:NIfTot,n_excits), source = temp_list(:,1:n_excits))

        call sort(det_list, ilut_lt, ilut_gt)

    end subroutine gen_all_singles_rs_hub_default

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
        integer(n_int), allocatable :: temp_dets(:,:)
        integer :: save_states

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

        if (tHPHF) then 
            ! for hphfs we need to purify the hilbert space! 
            save_states = n_states
            allocate(temp_dets(0:NIfTot,n_states), source = state_list_ilut)

            call spin_purify(save_states, temp_dets, n_states, state_list_ilut)

        end if

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

    subroutine pick_spin_opp_holes(ilutI, orbs, p_orb) 
        ! routine to pick two spin-opposite unoccupied orbitals 
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: orbs(2) 
        real(dp), intent(out) :: p_orb 
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "pick_spin_opp_holes"
#endif

        integer, parameter :: max_trials = 100
        integer :: cnt 

        cnt = 0
        orbs = 0
        p_orb = 0.0_dp 

        do while (cnt < max_trials)
            orbs(1) = 1 + int(genrand_real2_dsfmt() * nBasis) 

            if (IsOcc(ilutI,orbs(1))) cycle 

            do 
                orbs(2) = 1 + int(genrand_real2_dsfmt() * nBasis) 

                if (IsOcc(ilutI,orbs(2))) cycle 

                if (orbs(1) /= orbs(2)) exit 

            end do

            if (.not. same_spin(orbs(1),orbs(2))) exit 

        end do

        orbs = [minval(orbs),maxval(orbs)]

        p_orb = 1.0_dp / real((nBasis/2 - nOccAlpha)*(nBasis/2 - nOccBeta),dp)

    end subroutine pick_spin_opp_holes

    subroutine pick_random_hole(ilutI, orb, p_orb, spin) 
        ! routine to pick an unoccupied spin-orbital (orb) with probability 
        ! p_orb
        ! with an additional spin input we can restrict ourself to a 
        ! specific parallel spin! 
        ! spin = 1 -> beta orbital to be picked 
        ! spin = 0 -> alpha orbital to be picked
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: orb 
        real(dp), intent(out) :: p_orb
        integer, intent(in), optional :: spin 
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "pick_random_hole" 
#endif
        integer, parameter :: max_trials = 100
        integer :: cnt, i

        orb = 0
        cnt = 0
        p_orb = 0.0_dp

        if (present(spin)) then 
            ASSERT(spin == 0 .or. spin == 1)
          
            do while (cnt < max_trials)
                ! create a random spin-orbital of parallel spin
                i = 2*(1 + int(genrand_real2_dsfmt() * nBasis/2)) - spin

                if (IsNotOcc(ilutI,i)) then 
                    orb = i
                    if (spin == 0) then 
                        p_orb = 1.0_dp / real(nBasis/2 - nOccAlpha,dp)
                    else if (spin == 1) then 
                        p_orb = 1.0_dp / real(nBasis/2 - nOccBeta,dp)
                    end if
                    return
                end if
            end do

        else

            do while (cnt < max_trials)
                cnt = cnt + 1
                i = 1 + int(genrand_real2_dsfmt() * nBasis)

                if (IsNotOcc(ilutI, i)) then 
                    orb = i 
                    p_orb = 1.0_dp / real(nBasis - nel, dp)
                    return 
                end if
            end do
        end if

    end subroutine pick_random_hole

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

        integer(n_int) :: tmp_ilut(0:niftot)
        integer :: i, j, k, l, spin, ind
        

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

        if (tHPHF) then 
            ! i have to get the relevant HPHF determinant
            call EncodeBitDet(neel_state, tmp_ilut)
            tmp_ilut = return_hphf_sym_det(tmp_ilut)
            call decode_bit_det(neel_state, tmp_ilut)
        end if

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

    subroutine gen_all_excits_k_space_hubbard(nI, n_excits, det_list)!, sign_list) 

        integer, intent(in) :: nI(nel)
        integer, intent(out) :: n_excits 
        integer(n_int), intent(out), allocatable :: det_list(:,:)
        character(*), parameter :: this_routine = "gen_all_excits_k_space_hubbard"

        integer(n_int), allocatable :: triple_dets(:,:), temp_dets(:,:)
        integer :: n_triples, save_excits
        real(dp), allocatable :: sign_list(:)

        call gen_all_doubles_k_space(nI, n_excits, det_list, sign_list)

        if (t_trans_corr_2body) then 
            save_excits = n_excits
            ! also account for triple excitations
            call gen_all_triples_k_space(nI, n_triples, triple_dets)

            n_excits = n_excits + n_triples

            allocate(temp_dets(0:niftot, save_excits), source = det_list(:,1:save_excits))

            deallocate(det_list) 

            allocate(det_list(0:niftot,n_excits)) 

            det_list(:,1:save_excits) = temp_dets 

            det_list(:,save_excits+1:n_excits) = triple_dets

        end if

        call sort(det_list, ilut_lt, ilut_gt)

        ! if we have HPHF turned on we want to "spin-purify" the excitation 
        ! list 
        if (tHPHF) then 
            save_excits = n_excits

            if (allocated(temp_dets)) deallocate(temp_dets)

            allocate(temp_dets(0:NIfTot, n_excits), source = det_list)

            call spin_purify(save_excits, temp_dets, n_excits, det_list)

        end if

    end subroutine gen_all_excits_k_space_hubbard

    subroutine create_hilbert_space_kspace(n_alpha, n_beta, n_orbs, nI, & 
            n_states, state_list_ni, state_list_ilut) 
        ! new functionality to create all states of the hilber space in the 
        ! k-space hubbard model. reuse stuff from the real-space and apply 
        ! the momentum conservation criteria! 
        integer, intent(in) :: n_alpha, n_beta, n_orbs, nI(nel)
        integer, intent(out) :: n_states
        integer, intent(out), allocatable :: state_list_ni(:,:) 
        integer(n_int), intent(out), allocatable :: state_list_ilut(:,:) 
#ifdef __DEBUG
        character(*), parameter :: this_routine = "create_hilbert_space_kspace"
#endif

        integer(n_int), allocatable :: all_dets(:), temp_list_ilut(:,:) 
        integer :: i, j, nJ(nel), n_allowed, k_vec_in(3), k_vec(3)
        type(Symmetry) :: sym_prod, sym_in
        integer, allocatable :: temp_list_ni(:,:) 
        integer(n_int) :: temp_ilut(0:niftot)
        integer :: save_states

        ! this creates all possible basis states:
        all_dets = create_all_dets(n_orbs, n_alpha, n_beta) 

        sym_in = G1(nI(1))%Sym
        k_vec_in = G1(nI(1))%k
        do i = 2, nel 
            sym_in = SymTable(sym_in%s, G1(nI(i))%Sym%s)
            ! also add the k-vecs in a way so they do not leave the 
            ! first BZ 
            k_vec_in = lat%add_k_vec(k_vec_in, G1(nI(i))%k)

            ! a sanity check to see if it remained.. can be removed later! 
!             k_vec_in = lat%map_k_vec(k_vec_in + G1(nI(i))%k)
            
        end do

        allocate(temp_list_ilut(0:niftot, size(all_dets)))
        allocate(temp_list_ni(nel, size(all_dets)))

        temp_list_ilut = 0_n_int 
        temp_list_ni = 0
         
        ! now i want to loop over it to check the total momentum 
        ! how do i most efficiently determine the symmetry? 
        ! since adding the k-vector is not an option or? how far can the 
        ! momentum then go in terms of Brillouin zones? can i efficiently 
        ! map it back? or i just use the symmetry table, which has to be 
        ! setup beforehand! 
        n_allowed = 0

        do i = 1, size(all_dets) 
            temp_ilut = all_dets(i)
            call decode_bit_det(nJ, temp_ilut)
            sym_prod = G1(nJ(1))%sym
            k_vec = G1(nJ(1))%k
            do j = 2, nel 
                ! and i have to decode into the nI representation.. 
                ! oh god this is inefficient.. 
                sym_prod = symtable(sym_prod%s, G1(nJ(j))%Sym%s)
                k_vec = lat%add_k_vec(k_vec, G1(nJ(j))%k)

!                 k_vec = lat%map_k_vec(k_vec + G1(nJ(j))%k)
            end do
            if (sym_prod%s == sym_in%s) then 
                ! if the symmetry fits: 
                ASSERT(all(k_vec == k_vec_in))
                n_allowed = n_allowed + 1
                temp_list_ilut(:,n_allowed) = temp_ilut
                temp_list_ni(:,n_allowed) = nJ
            else 
                ASSERT(.not. all(k_vec == k_vec_in))
            end if
        end do

        n_states = n_allowed

        allocate(state_list_ni(nel,n_allowed), source = temp_list_ni(:,1:n_allowed))
        allocate(state_list_ilut(0:niftot,n_allowed), source = temp_list_ilut(:,1:n_allowed))

        if (tHPHF) then 
            ! for hphfs we need to purify the hilbert space! 
            save_states = n_states
            if (allocated(temp_list_ilut)) deallocate(temp_list_ilut)

            allocate(temp_list_ilut(0:NIfTot,n_states), source = state_list_ilut)

            call spin_purify(save_states, temp_list_ilut, n_states, state_list_ilut)

        end if

    end subroutine create_hilbert_space_kspace

    function get_orb_from_kpoints_three(src, orb_a, orb_b) result(orb_c) 
        ! momentum conservation for three-body terms 
        integer, intent(in) :: src(3), orb_a, orb_b
        integer :: orb_c 
#ifdef __DEBUG
        character(*), parameter :: this_routine = "get_orb_from_kpoints_three"
#endif
        integer :: sum_ms, kc(3), spin_c, spin_ab
        type(symmetry) :: kc_sym

        ! implement that generally for also all-spin parallel excitation, which 
        ! might be necessary in the future.. 
        sum_ms = sum(get_spin_pn(src))

        ASSERT(sum_ms == -3 .or. sum_ms == -1 .or. sum_ms == 1 .or. sum_ms == 3)

        ! momentum conservation: ka + kb + kc = ki + kj + kk
        ! to this better with the symmetry table or a new addition 
        ! functionality for k-vectors, to never leave the first BZ
!         kc = G1(src(1))%k + G1(src(2))%k + G1(src(3))%k - G1(orb_a)%k - G1(orb_b)%k 
        ! do something like: 
        ! and write two routines, interfaced, one take a vector of 
        ! k-values the other takes the encoded symmetry symbol
        ! and setup a symmetry table and inverse symmetry list like in the 
        ! sym mod
!         kc = lat%add_k_vec(lat%add_k_vec(G1(src(1))%k, G1(src(2))%k), G1(src(3))%k)
        kc_sym = SymTable(G1(src(1))%sym%s, SymTable(G1(src(2))%sym%s, G1(src(3))%sym%s)%s)
!         kc = lat%subtract_k_vec(kc, lat%add_k_vec(G1(orb_a)%k,G1(orb_b)%k))
        kc_sym = SymTable(kc_sym%s, SymConjTab(SymTable(G1(orb_a)%sym%s,G1(orb_b)%sym%s)%s))

        kc = lat%sym_to_k(kc_sym%s,:)

        !do it over the mult table now!!

        ! perdiodic BC: 
        if (tHub) then
            if (t_k_space_hubbard) then 
                ! use the new lattice implementation
                ! with the above correct addition it should still be in the 
                ! first BZ.. 
!                 kc = lat%map_k_vec(kc)
            else
                call mompbcsym(kc, nBasisMax) 
            end if
        end if

        ! and now we want the spin: 
        if (sum_ms == -3) then 
            ! assert we can still reach the desired spin:
            ASSERT(get_spin_pn(orb_a) + get_spin_pn(orb_b) == -2)

            spin_c = 1 
        else if (sum_ms == 3) then 
            ASSERT(get_spin_pn(orb_a) + get_spin_pn(orb_b) == 2) 

            spin_c = 2

        else 
            spin_ab = get_ispn([orb_a,orb_b])

            ! we know we are not totally parallel so if: 
            if (spin_ab == 1) then 
                ! if we have already two beta, we need alpha
                spin_c = 2
            else if (spin_ab == 3) then 
                ! if we have two alpha we need a beta 
                spin_c = 1 

            else 
                ! in this case we need a spin to reach the final ms
                if (sum_ms == -1) then 
                    spin_c = 1
                else 
                    spin_c = 2
                end if
            end if
        end if

        if (t_k_space_hubbard) then 
            orb_c = lat%get_orb_from_k_vec(kc, spin_c)
        else
            orb_c = KPointToBasisFn(kc(1),kc(2),kc(3),spin_c) 
        end if

    end function get_orb_from_kpoints_three

    subroutine pick_three_opp_elecs(nI, elecs, p_elec, opt_sum_ms) 
        integer, intent(in) :: nI(nel)
        integer, intent(out) :: elecs(3)
        real(dp), intent(out) :: p_elec 
        integer, intent(out), optional :: opt_sum_ms
#ifdef __DEBUG
        character(*), parameter :: this_routine = "pick_three_opp_elecs"
#endif
        integer :: sum_ms

        ! this can be done way more effective than this simple implementation:
        first: do 
            elecs(1) = 1 + int(genrand_real2_dsfmt() * nel) 

            do 
                elecs(2) = 1 + int(genrand_real2_dsfmt() * nel) 

                if (elecs(1) /= elecs(2)) exit 
            end do

            if (get_ispn(nI(elecs(1:2))) /= 2) then 
                ! if the already picked electrons have same spin
                ! take opposite spin
                do 
                    elecs(3) = 1 + int(genrand_real2_dsfmt() * nel) 

                    if (elecs(1) /= elecs(3) .and. elecs(2) /= elecs(3) .and. &
                        .not. same_spin(nI(elecs(1)),nI(elecs(3)))) then 
                        exit first
                    end if
                end do
            else 
                ! here we can pick any spin electron 

                do 
                    elecs(3) = 1 + int(genrand_real2_dsfmt() * nel) 

                    if (elecs(1) /= elecs(3) .and. elecs(2) /= elecs(3)) then 
                        exit first
                    end if 
                end do
            end if
        end do first
                        
        ! i have to be careful how i could have gotten the electrons.. 
        ! because the order does not matter and there is no restriction set 
        ! output an ordered form
        elecs = sort_unique(elecs)

        ! the first two electrons are picked totally randon, independent of 
        ! spin just with the restriction i /= j 
        ! so p(i) = 1/nel, p(j|i) = 1/(nel-1) 
        ! the third electron has a spin restriction.. 
        ! if spin(i) = spin(j) then spin(k) must be -spin(i) giving 
        ! p(k|ij) = 1/n_opp_spin 
        ! if spin(i) /= spin(j) then k is freely except /= i,j 
        ! and since the order is not important the total probability is 
        ! p(i) * p(i|j) * [p(k|ij) + p(i|jk) + p(j|ik)] 
        ! which translates to 
        ! 1/nel * 1/(nel-1) * [1/n_opp + 2/(nel-2)] 
        ! so depending on the total mz, we get 

        sum_ms = sum(get_spin_pn(nI(elecs)))
        ASSERT(sum_ms == 1 .or. sum_ms == -1)
        if (sum_ms == 1) then 
            ! then we have 2 alpha and 2 beta electrons 
            p_elec = 2.0_dp/real(nel*(nel-1),dp) * &
                (1.0_dp/real(nOccBeta,dp) + 2.0_dp/real(nel-2,dp))
        else if (sum_ms == -1) then
            ! then we have 2 beta electrons and 1 alpha 
            p_elec = 2.0_dp/real(nel*(nel-1),dp) * & 
                (1.0_dp/real(nOccAlpha,dp) + 2.0_dp/real(nel-2,dp))
        end if

        ! adjust the edge cases 
!         if (nOccAlpha == 1 .or. nOccBeta == 1) p_elec = 2.0_dp * p_elec

        if (present(opt_sum_ms)) opt_sum_ms = sum_ms

    end subroutine pick_three_opp_elecs

    subroutine pick_spin_par_elecs(nI, elecs, p_elec, opt_ispn)
        integer, intent(in) :: nI(nel) 
        integer, intent(out) :: elecs(2)
        real(dp), intent(out) :: p_elec 
        integer, intent(out), optional :: opt_ispn
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "pick_spin_par_elecs"
#endif 
        integer :: ispn

        do 
            elecs(1) = 1 + int(genrand_real2_dsfmt() * nel)

            do 
                elecs(2) = 1 + int(genrand_real2_dsfmt() * nel)
                
                if (elecs(1) /= elecs(2)) exit 

            end do
            ispn = get_ispn(nI(elecs))
            if (ispn /= 2) exit 

        end do

        ASSERT(ispn == 1 .or. ispn == 3) 

        ! and always output an ordered form 
        elecs = [minval(elecs),maxval(elecs)]

        ! i have to rework the probability here.. 
        if (ispn == 1) then 
            p_elec = 1.0_dp / real(nOccBeta*(nOccBeta-1),dp)
        else if (ispn == 3) then 
            p_elec = 1.0_dp / real(nOccAlpha*(nOccAlpha-1),dp)
        end if

        ! in the special case of all spin-parallel:
        if (nOccBeta == 0 .or. nOccAlpha == 0) p_elec = 2.0_dp * p_elec

        if (present(opt_ispn)) opt_ispn = ispn

    end subroutine pick_spin_par_elecs

    function find_minority_spin(spins) result(opp) 
        ! for now, given 3 spin orbitals, where 2 share a spin, this function
        ! returns the opposite spin
        integer, intent(in) :: spins(:) 
        integer :: opp
#ifdef __DEBUG
        character(*), parameter :: this_routine = "find_minority_spin" 
#endif 
        integer :: i

        ! for now, do it only for 3 inputted spins:
        ASSERT(size(spins) == 3)
        ASSERT(sum(get_spin_pn(spins)) == -1 .or. sum(get_spin_pn(spins)) == 1)

        opp = -1

        if (sum(get_spin_pn(spins)) == -1) then 
            ! then we have 2 beta and 1 alpha 
            do i = 1, size(spins)
                if (is_alpha(spins(i))) opp = spins(i)
            end do
        else
            ! we have 2 alpha and 1 beta 
            do i = 1, size(spins)
                if (is_beta(spins(i))) opp = spins(i)
            end do
        end if
    end function find_minority_spin

    subroutine gen_all_triples_k_space(nI, n_excits, det_list) 
        integer, intent(in) :: nI(nel) 
        integer, intent(out) :: n_excits 
        integer(n_int), intent(out), allocatable :: det_list(:,:)
        character(*), parameter :: this_routine = "gen_all_triples_k_space"

        integer :: i,j,k,a,b,c,spin,ind,src(3),ex(2,3),n_bound,pos
        integer(n_int) :: ilut(0:niftot), ilutJ(0:NIfTot) 
        integer(n_int), allocatable :: temp_list(:,:)
        real(dp) :: elem

        ! write a routine to create all the triple excitations for the 
        ! transcorrelated k-space hubbard model. for a given determinant 
        call EncodeBitDet(nI, ilut)

        ! do it in the most naive way for now.. 
        ! i need to loop over all different electron combinations, with the 
        ! restriction that they are not all parallel spin (due to the form 
        ! of the transcorrelated hamiltonian) 
        n_excits = 1

        ! todo: an estimate for the upper bound of number of triple excitations.. 
        ! this gets too high for big lattices.. 

!         n_bound = nel*(nel-1)*(nel-2) * (nbasis - nel)*(nbasis - nel -1)*(nbasis-nel-2)
        ! i think a more correct estimat is:
        n_bound = int(nel*(nel-1)*(nel-2) * (nbasis - nel)*(nbasis - nel - 1)/8)

        allocate(temp_list(0:niftot,n_bound))
        temp_list = 0_n_int

        do i = 1, nel - 2
            do j = i + 1, nel - 1
                do k = j + 1, nel 
                    ! assert that the they are not all spin-parallel: 
                    src = nI([i,j,k])
                    spin = sum(get_spin_pn(src))
                    ex(1,:) = src

                    if (spin == -1 .or. spin == 1) then 
                        do a = 1, nbasis - 1
                            if (IsNotOcc(ilut,a)) then 
                                do b = a + 1, nbasis 
                                    if (IsNotOcc(ilut,b)) then 

                                        ! this gives me the correct spin already 
                                        ! to have the same as the electrons!
                                        c = get_orb_from_kpoints_three(src, a, b) 

                                        ! i also have to check if no orbital is 
                                        ! picked twice
                                        if (IsNotOcc(ilut, c) .and. c/= a .and. c /=b) then 
                                            ! this should be a valid excitation or?

                                            ex(2,:) = [a,b,c] 
                                            ! should i check the matrix element too?
                                            ! and also be sure that the momentum fits
!                                             if (abs(get_3_body_helement_ks_hub(nI,ex,.false.)) > EPS) then 
                                            elem = abs(get_helement_lattice(nI, 3, ex, .false.))
                                            if (elem > EPS) then
                                                ilutJ = make_ilutJ(ilut, ex,3)

                                                ! and to be save, search if we have 
                                                ! this excitation already.. 
                                                pos = binary_search(temp_list(:,1:n_excits), ilutJ, nifd+1)

                                                if (pos < 0) then 
                                                    ! new excitation
                                                    temp_list(:,n_excits) = ilutJ
                                                    n_excits = n_excits + 1 
                                                end if
                                            end if
                                        end if
                                    end if
                                end do
                            end if
                        end do
                    end if
                end do
            end do
        end do

        n_excits = n_excits - 1
        allocate(det_list(0:NIfTot,n_excits), source=temp_list(:,1:n_excits))

    end subroutine gen_all_triples_k_space

    subroutine gen_all_doubles_k_space(nI, n_excits, det_list, sign_list)
        integer, intent(in) :: ni(nel) 
        integer, intent(out) :: n_excits
        integer(n_int), intent(out), allocatable :: det_list(:,:) 
        real(dp), intent(out), allocatable, optional :: sign_list(:)
        character(*), parameter :: this_routine = "gen_all_doubles_k_space"

        integer(n_int) :: ilutJ(0:niftot), ilut(0:niftot)
        integer :: n_bound, i, j, a, b, src(2), ex(2,2), pos, n_par, n_opp, nj(nel)
        integer(n_int), allocatable :: temp_list(:,:) 
        real(dp) :: elem
        logical :: t_sign, tpar
        real(dp), allocatable :: temp_sign(:)

        ! fuck it! these annoying old routines break my balls! 
        ! just write a new one to test my excitations with! 
        call EncodeBitDet(nI, ilut) 
        
        n_excits = 1 

!         n_bound = nel*(nel-1)*(nbasis-nel)*(nbasis-nel-1)
        ! i think a more correct estimat is:
        n_bound = max(int(nel*(nel-1)*(nBasis - nel)/4), 10)

        allocate(temp_list(0:niftot,n_bound))

        n_par = 0
        n_opp = 0

        temp_list = 0_n_int 

        if (present(sign_list)) then 
            t_sign = .true. 
            allocate(temp_sign(n_bound), source = 0.0_dp)
        else 
            t_sign = .false.
        end if

        do i = 1, nel -1 
            do j = i + 1, nel 
                src = nI([i,j])
                
                ex(1,:) = src 

                ! if we do not have transcorrelation cyclce for same spins
                if (same_spin(src(1),src(2)) .and. .not. t_trans_corr_2body) cycle 

                do a = 1, nbasis
                    if (IsNotOcc(ilut,a)) then 
                        b = get_orb_from_kpoints(src(1),src(2), a) 

                        if (IsNotOcc(ilut,b) .and. a /= b) then 

                            ! do i need to sort ex?? try..
                            ex(2,:) = [a,b] 

!                             if (abs(get_offdiag_helement_k_sp_hub(nI, ex, .false.)) > EPS) then 
                            ! use the get_helement_lattice to avoid circular 
                            ! dependencies
                            call make_double(nI,nJ,i,j,a,b,ex,tpar)
                            elem = get_helement_lattice(nI,nJ)
                            !elem = get_helement_lattice(nI, 2, ex, .false.)
                            if (abs(elem) > EPS) then

                                ilutJ = make_ilutJ(ilut, ex, 2) 

                                pos = binary_search(temp_list(:,1:n_excits), ilutJ,nifd+1)

                                if (pos < 0) then 

                                    
                                    if (t_sign) then 
                                         temp_sign(n_excits) = sign(1.0_dp, elem)
                                         
                                    end if
                                    temp_list(:,n_excits) = ilutJ
                                    !call sort(temp_list(:,1:n_excits),ilut_lt,ilut_gt)
                                    call sort(temp_list(:,1:n_excits),temp_sign(1:n_excits))

                                    n_excits = n_excits + 1
                                    ! damn.. i have to sort everytime i guess..

                                    ! to be sure also count seperately if it is 
                                    ! a parallel or opposite excitation 
                                    if (same_spin(src(1),src(2))) then 
                                        ASSERT(same_spin(src(1),a))
                                        ASSERT(same_spin(a,b))
                                        n_par = n_par + 1 
                                    else 
                                        ASSERT(.not. same_spin(a,b))
                                        n_opp = n_opp + 1 
                                    end if
                                end if
                            end if
                        end if
                    end if
                end do
            end do
        end do

        n_excits = n_excits - 1
        allocate(det_list(0:NIfTot,n_excits), source = temp_list(:,1:n_excits))

        if (t_sign) then 

            allocate(sign_list(n_excits), source = temp_sign(1:n_excits))
	    !print *, "before:"
	    !do i = 1, n_excits
		!print *, det_list(:,i), sign_list(i)
 	    !end do
	    !print *, "after:"
            call sort(det_list, sign_list)!, ilut_lt, ilut_gt)
	    !do i = 1, n_excits
		!print *, det_list(:,i), sign_list(i)
 	    !end do
        else 
            call sort(det_list, ilut_lt, ilut_gt)
        end if

    end subroutine gen_all_doubles_k_space

    function get_orb_from_kpoints(orbi, orbj, orba) result(orbb)
        ! write a cleaner implementation of this multiple used 
        ! functionality.. because kPointToBasisFn to basisfunction 
        ! promises more than it actually odes 
        integer, intent(in) :: orbi, orbj, orba
        integer :: orbb

        integer :: ki(3), kj(3), ka(3), kb(3), ispn, spnb


        ki = G1(orbi)%k
        kj = G1(orbj)%k

        ! Given A, calculate B in the same way as before
        ka = G1(orba)%k
        kb = ki + kj - ka

        ispn = get_ispn([orbi,orbj])

        if (iSpn == 1) then
            spnb = 1
        elseif (iSpn == 2) then
            ! w.d: bug found by peter jeszenski and confirmed by 
            ! simon! 
            ! i do not think this is so much more efficient than an 
            ! additional if-statement which is way more clear!
            ! messed up alpa and beta spin here.. 
!             spnb = (-G1(orba)%Ms + 1)/2 + 1
    !       spnb = (G1(orba)%Ms)/2 + 1
            if (is_beta(orba)) then 
                spnb = 2
            else 
                spnb = 1
            end if
        elseif(iSpn == 3) then
            spnb = 2
        end if

        ! damn.. for some reason this is different treated in the 
        ! hubbard and UEG case..
        ! i turn off the thub flag in my new implementation, so this is 
        ! never actually reached below.. i should change that 

        if (t_k_space_hubbard) then 
            orbb = lat%get_orb_from_k_vec(kb, spnb)
            return
        end if

        ! this now distinguishes between UEG and hubbard models.
        if (tHub) then
            call mompbcsym(kb, nbasismax)
        end if

        orbb = KPointToBasisFn(kb(1), kb(2), kb(3), spnb)

    end function get_orb_from_kpoints

    function get_ispn(src) result(ispn)
        ! it is annoying to write this over and over again.. 
        integer, intent(in) :: src(2) 
        integer :: ispn 

        if (is_beta(src(1)) .eqv. is_beta(src(2))) then 
            if (is_beta(src(1))) then 
                ispn = 1
            else 
                ispn = 3
            end if
        else 
            ispn = 2 
        end if
    end function get_ispn

    function make_ilutJ(ilutI, ex, ic) result(ilutJ)
        ! function similar to make_single and make_double to create the 
        ! accompaning ilut form. 
        integer(n_int), intent(in) :: ilutI(0:niftot) 
        integer, intent(in) :: ic
        integer, intent(in) :: ex(2,ic)
        integer(n_int) :: ilutJ(0:niftot) 

#ifdef __DEBUG
        character(*), parameter :: this_routine = "make_ilutJ"
#endif

        integer :: ij(ic), ab(ic), i

#ifdef __DEBUG
        ASSERT(ic == 1 .or. ic == 2 .or. ic == 3)
        ! should this every be called with 0 orbitals.. i guess no.. 
        do i = 1, ic 
            ASSERT(ex(1,ic) > 0) 
            ASSERT(ex(2,ic) > 0)
            ASSERT(ex(1,ic) <= nbasis)
            ASSERT(ex(2,ic) <= nbasis)
        end do
#endif

        ilutJ = ilutI 
        

        ij = get_src(ex)
        ab = get_tgt(ex)

        do i = 1, ic 
            clr_orb(ilutJ, ij(i))
            set_orb(ilutJ, ab(i))
        end do

    end function make_ilutJ

    function find_elec_in_ni(nI, orb) result(elec)
        ! routine to find the number of the elctron in spin-orbital orb
        integer, intent(in) :: nI(nel), orb
        integer :: elec 
#ifdef __DEBUG
        character(*), parameter :: this_routine = "find_elec_in_ni"
#endif

        ASSERT(orb > 0)
        ASSERT(orb <= nbasis) 

        ! can i just reuse 
        elec = binary_search_first_ge(nI, orb)

        ! it already make the fail case.. or?
        ! and then make a fail-case: 
        if (elec == -1) return 

        if (nI(elec) /= orb) then 
            ! it is actually not found? 
            elec = -1 
        end if

    end function find_elec_in_ni

    function get_occ_neighbors(ilut, orb) result(occ_neighbors) 
        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer, intent(in) :: orb 
        real(dp) :: occ_neighbors
#ifdef __DEBUG
        character(*), parameter :: this_routine = "get_occ_neighbors" 
#endif
        integer, allocatable :: neighbors(:)
        integer :: i

        ASSERT(associated(lat)) 

        ! orb is given as a spatial orbital! 

        neighbors = lat%get_neighbors(orb) 

        occ_neighbors = 0.0_dp
        do i = 1, size(neighbors)
            ! check both spinorbitals
            if (IsOcc(ilut, 2*neighbors(i)))   occ_neighbors = occ_neighbors + 1.0_dp
            if (IsOcc(ilut, 2*neighbors(i)-1)) occ_neighbors = occ_neighbors + 1.0_dp
        end do

    end function get_occ_neighbors

    function get_spin_density_neighbors(ilut, orb) result(spin_dens_neighbors)
        ! function to get the spin-density of the neighboring orbitals 
        ! n_{i,\up} - n_{i,\down} 
        integer(n_int), intent(in) :: ilut(0:NIfTot) 
        integer, intent(in) :: orb
        real(dp) :: spin_dens_neighbors
#ifdef __DEBUG
        character(*), parameter :: this_routine = "get_spin_density_neighbors"
#endif 
        integer, allocatable :: neighbors(:) 
        integer :: i 

        ASSERT(associated(lat)) 

        spin_dens_neighbors = 0.0_dp 

        ! orb is given in spatial orbials 
        neighbors = lat%get_neighbors(orb) 
        do i = 1, size(neighbors)
            ! what do we degine as up spin?? have to be sure here
            ! lets take alpha
            if (IsOcc(ilut,2*neighbors(i)-1)) spin_dens_neighbors = spin_dens_neighbors - 1.0_dp 
            if (IsOcc(ilut,2*neighbors(i)))   spin_dens_neighbors = spin_dens_neighbors + 1.0_dp
        end do

        spin_dens_neighbors = spin_dens_neighbors / 2.0_dp

    end function get_spin_density_neighbors


end module lattice_models_utils
