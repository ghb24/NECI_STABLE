#include "macros.h"
! small utitlities to analyse the point group symmetry of the 
! FCIQMC wavefunction 

module analyse_wf_symmetry

    use SystemData, only: nBasis, g1, arr, brr, nel

    use Parallel_neci, only: iProcIndex, root, MPI_2DOUBLE_PRECISION, &
                             MPI_MAXLOC, MPIAllReduceDatatype, MPIBCast, &
                             MPISum

    use fcimcdata, only: CurrentDets, TotWalkers

    use semi_stoch_procs, only: return_most_populated_states

    use lattice_mod, only: lat

    use bit_rep_data, only: niftot, nifd

    use constants, only: n_int, dp, pi, lenof_sign

    use util_mod, only: binary_search, binary_search_int

    use bit_reps, only: extract_sign, encode_sign, decode_bit_det, &
                        encode_bit_rep

    use DetBitOps, only: EncodeBitDet, ilut_lt, ilut_gt, DetBitEq

    use sort_mod, only: sort

    use unit_test_helpers, only: print_matrix

    implicit none

    logical :: t_symmetry_analysis = .false.
    logical :: t_symmetry_rotation = .false. 
    real(dp) :: symmetry_rotation_angle = 0.0_dp 

    logical :: t_symmetry_mirror = .false. 
    character(1) :: symmertry_mirror_axis = '0'

    logical :: t_symmetry_inversion = .false.

    logical :: t_read_symmetry_states = .false. 
    integer :: n_symmetry_states = 0
    logical :: t_pop_symmetry_states = .false. 

    integer, allocatable :: symmetry_states(:,:)
    real(dp), allocatable :: symmetry_weights(:)
    integer(n_int), allocatable :: symmetry_states_ilut(:,:)

    interface rotate
        module procedure rotate_orb
        module procedure rotate_vec
    end interface rotate

contains

    subroutine analyze_wavefunction_symmetry()
        ! driver routine to analyze the symmetry of a part of the 
        ! wavefunction. 
        integer :: sym_orbs(nBasis/2), orig_orbs(nBasis/2)
        integer :: i, transformed_states(nel, n_symmetry_states)
        real(dp) :: transformed_weights(n_symmetry_states)
        integer(n_int) :: transformed_states_ilut(0:niftot,n_symmetry_states)

        ! i have the original list of the orbitals and the 
        ! corresponding k- or r-vectors. 

        call init_symmetry_states()

        if (iProcIndex == root) then
            print *, "Analyze the symmetry of the FCIQMC wavefunction: "

            if (t_symmetry_rotation) then 
                print *, "applying rotation of ", symmetry_rotation_angle, " degrees"
            end if

            if (t_symmetry_mirror) then 
                print *, "mirror wf. along ", symmertry_mirror_axis, "-axis"
            end if

            if (t_symmetry_inversion) then
                print *, "apply inversion symmetry"
            end if

            if (t_read_symmetry_states) then
                print *, "using the ", n_symmetry_states, " read-in states: "
                do i = 1, n_symmetry_states
                    print *, symmetry_states(:,i)
                end do
            end if

            if (t_pop_symmetry_states) then 
                print *, "using the ", n_symmetry_states,  " highest occupied states"
                do i = 1, n_symmetry_states
                    print *, symmetry_states(:,i)
                end do
            end if

            ! first check if everything is setup correctly.. 
            print *, "brr:" 
            do i = 1, nBasis
                print *, brr(i)
            end do

            print *, "i, lat%ind, lat%get_sym, lat%get_k_vec"
            do i = 1, nBasis/2
                print *, i, lat%get_site_index(i), lat%get_sym(i), lat%get_k_vec(i)
            end do

            call WRITEBASIS(6,g1,nBasis,arr,brr)

            ! i need to apply the chosen symmetry to the vectors and determine, 
            ! which orbital maps into which 
            ! the flag if we actually apply is within the functions.
            ! ROTATION:
            orig_orbs = get_spatial(brr(1:nBasis:2))
            sym_orbs = apply_rotation(orig_orbs)

            
            print *, "orig orbs -> sym_orbs: "
            do i = 1, nBasis/2
                print *, orig_orbs(i), " -> ", sym_orbs(i)
            end do

            ! i can ofc apply multiple symmetries
            ! MIRROR
    !         sym_orbs = apply_mirror(sym_orbs)

            ! INVERSION
    !         sym_orbs = apply_inversion(sym_orbs)

            ! now i have the mapping between the original and the 
            ! symmetry transformed orbitals

            call transform_states(orig_orbs, sym_orbs, transformed_states, &
                transformed_weights, transformed_states_ilut)

            call compare_states(symmetry_states_ilut, transformed_states_ilut)

        end if

    end subroutine analyze_wavefunction_symmetry

    subroutine compare_states(orig_states, trans_states)
        ! compare the original and the symmetry transformed states 
        integer(n_int), intent(in) :: orig_states(0:niftot, n_symmetry_states), &
                                      trans_states(0:niftot, n_symmetry_states)

        integer(n_int) :: null_int(0:niftot), ilutI(0:niftot), ilutJ(0:niftot)
        integer :: i, j, k, l
        ! print the original and the transformed states next to each other 
        ! and missing determinants in the repective lists are indicated 
        i = 1 
        j = 1

        null_int = 0_n_int

        ! or first create a list and an integer indicator, where the 
        ! determinant is from.. 
        do k = 1, 2*n_symmetry_states
            ilutI = orig_states(:,i)
            ilutJ = trans_states(:,j)

            if (DetBitEq(ilutI,ilutJ)) then 
                ! if both are equal we can move on 
                call print_2_states(ilutI,ilutJ)
                i = i + 1
                j = j + 1

            else if (ilut_lt(ilutI,ilutJ)) then 
                ! if I is less than J, we have to increase I and only print I
                call print_2_states(ilutI, null_int)
                i = i + 1

            else 
                ! if J is smaller we have to print J and increase J
                call print_2_states(null_int, ilutJ)
                j = j + 1

            end if
            
            ! and provide the correct exci conditions
            if (i == j .and. i > n_symmetry_states) exit

            if (i > n_symmetry_states .and. j <= n_symmetry_states) then 
                ! if i is already above the list then the rest of the 
                ! entries are from list J
                do l = j, n_symmetry_states
                    ilutJ = trans_states(:,l)
                    call print_2_states(null_int, ilutJ)
                end do
                exit
            end if
            ! and if j is alreay above, then the rest is from I 
            if (i <= n_symmetry_states .and. j > n_symmetry_states) then 
                do l = i, n_symmetry_states
                    ilutI = orig_states(:,l)
                    call print_2_states(ilutI, null_int)
                end do
                exit 
            end if
        end do

    end subroutine compare_states

    subroutine print_2_states(left, right, nunit)
        integer(n_int), intent(in) :: left(0:niftot), right(0:niftot)
        integer, intent(in), optional :: nunit

        integer :: iout
        real(dp) :: left_sign(lenof_sign), right_sign(lenof_sign)

        if (present(nUnit)) then 
            iout = nunit
        else 
            iout = 6
        end if

        call extract_sign(left, left_sign)
        call extract_sign(right, right_sign)
        ! print the left side
        if (all(left == 0_n_int)) then 
            call print_null_det(iout)
        else
            call writeDetBit(iout, left,.false.)
        end if
        write(iout, '(A1)', advance = 'no') "|"
        write(iout, '(f16.7)', advance = 'no') left_sign(1)
        write(iout, '(A1)', advance = 'no') "|"
        write(iout, '(f16.7)', advance = 'no') right_sign(1)
        write(iout, '(A1)', advance = 'no') "|"


        ! then print right side 
        if (all(right == 0_n_int)) then 
            call print_null_det(iout,.true.)
        else
            call writeDetBit(iout, right,.true.)
        end if

    end subroutine print_2_states

    subroutine print_null_det(nunit, t_newline)
        integer, intent(in) :: nunit
        logical, intent(in), optional :: t_newline
        
        integer :: i 

        do i = 1, nBasis - 1
           write(nunit, "(A1)", advance = 'no') "-" 
       end do

       if (present(t_newline) .and. t_newline) then 
           write(nunit, "(A1)") "-"
       else
           write(nunit, "(A1)", advance = 'no') "-"
       end if

    end subroutine print_null_det

    subroutine transform_states(orig_orbs, transformed_orbs, &
            transformed_states, transformed_weights, transformed_states_ilut)
        integer, intent(in) :: orig_orbs(nBasis/2), transformed_orbs(nBasis/2)
        integer, intent(out) :: transformed_states(nel, n_symmetry_states)
        real(dp), intent(out) :: transformed_weights(n_symmetry_states)
        integer(n_int), intent(out) :: transformed_states_ilut(0:niftot,n_symmetry_states)
#ifdef __DEBUG
        character(*), parameter :: this_routine = "transform_states"
#endif
        integer :: i, n_phase, ind(n_symmetry_states)
        real(dp) :: tmp_sign(lenof_sign)
    
        tmp_sign = 0.0_dp

        ! do it plain ans stupid for now.. 
        do i = 1, n_symmetry_states

             call apply_transformation(symmetry_states(:,i), orig_orbs, &
                 transformed_orbs, transformed_states(:,i), n_phase)

             transformed_weights(i) = real(n_phase,dp) * symmetry_weights(i)

             tmp_sign(1) = transformed_weights(i)

             call EncodeBitDet(transformed_states(:,i), transformed_states_ilut(:,i))
             call encode_sign(transformed_states_ilut(:,i), tmp_sign)

        end do

        ! the original highest pop list is sorted by weight i guess.. 
        ! sort them by the integers in the ilut representation
        ind = [(i, i = 1, n_symmetry_states)]

        call sort(symmetry_states_ilut, ind)

        symmetry_states = symmetry_states(:,ind)
        symmetry_weights = symmetry_weights(ind)


        ind = [(i, i = 1, n_symmetry_states)]

        call sort(transformed_states_ilut, ind)

        transformed_states = transformed_states(:,ind)
        transformed_weights = transformed_weights(ind)

        
    end subroutine transform_states

    subroutine apply_transformation(nI, orig_orbs, trans_orbs, nJ, n_phase)
        ! apply the transformation encoded in the orig_orbs and trans_orbs 
        ! to nI to obtain nJ. nJ will not be sorted except
        ! the phase is present, where the fermionic phase from reshuffling the 
        ! orbitals into natural order is also outputted. 
        integer, intent(in) :: nI(nel), orig_orbs(nBasis/2), trans_orbs(nBasis/2)
        integer, intent(out) :: nJ(nel)
        integer, intent(out), optional :: n_phase
#ifdef __DEBUG
        character(*), parameter :: this_routine = "apply_transformation"
#endif
        integer :: i, pos

        nJ = 0

        do i = 1, nel 

            ! i cant use binary search, since orig orbs is not sorted! 
!             pos = binary_search_int(orig_orbs, get_spatial(nI(i)))
            pos = stupid_search(orig_orbs, get_spatial(nI(i)))

            ! it has to be found! 
            ASSERT(pos > 0)
            if (is_beta(nI(i))) then 
                nJ(i) = 2*trans_orbs(pos) - 1
            else 
                nJ(i) = 2*trans_orbs(pos)
            end if

        end do

        if (present(n_phase)) then 
            call sort(nJ, par = n_phase)
        end if 

    end subroutine apply_transformation

    function stupid_search(list, val) result(pos) 
        integer, intent(in) :: list(:), val
        integer :: pos 

        integer :: i

        pos = 0

        do i = lbound(list,1), ubound(list,1)
            if (val == list(i)) then 
                pos = i
                return
            end if
        end do

    end function stupid_search

    subroutine init_symmetry_states()
        ! routine to initialize the to be analysed states of the 
        ! wavefunction 
        integer :: i, n_states
        integer(n_int) :: ilut(0:niftot)
        integer(n_int), allocatable :: largest_walkers(:,:)
        real(dp) :: temp_sign(lenof_sign)

        if (t_read_symmetry_states) then 
            ! if we read them in, we have to find it in the walker-list 
            ! but only in the 1000 most populated states
            n_states = 1000

            do i = 1, n_symmetry_states
                call EncodeBitDet(symmetry_states(:,i), symmetry_states_ilut(:,i))
            end do

        else if (t_pop_symmetry_states) then 
            ! here we have to take the N most occupied from the 
            ! wavefunction 
            ! otherwise just find the N most populated states
            n_states = n_symmetry_states

        end if

        allocate(largest_walkers(0:niftot, n_states))
        largest_walkers = 0_n_int

        call get_highest_pop(n_states, largest_walkers)

        if (t_read_symmetry_states) then 
            call find_states_in_list(n_symmetry_states, symmetry_states, &
                largest_walkers, symmetry_weights)

        else if (t_pop_symmetry_states) then

            symmetry_states_ilut = largest_walkers(:,1:n_symmetry_states)
            
            do i = 1, n_symmetry_states

                call decode_bit_det(symmetry_states(:,i), largest_walkers(:,i))
                call extract_sign(largest_walkers(:,i), temp_sign)
                symmetry_weights(i) = temp_sign(1)

            end do
        end if

    end Subroutine init_symmetry_states

    subroutine find_states_in_list(n_states, nI_search, ilut_list, nI_weights)
        ! routine to find the nI in the ilut_list and assign the 
        ! corresponding weights. if not found the weights are 0
        integer, intent(in) :: n_states
        integer, intent(in) :: nI_search(nel,n_states)
        integer(n_int), intent(in) :: ilut_list(0:niftot,n_states)
        real(dp), intent(out) :: nI_weights(n_states)
#ifdef __DEBUG 
        character(*), PARAMETER :: this_routine = "find_states_in_list"
#endif 
        integer(n_int) :: ilut(0:niftot)
        integer :: i, pos
        real(dp) :: temp_sign(lenof_sign)

        nI_weights = 0.0_dp

        do i = 1, n_states 
            call EncodeBitDet(nI_search(:,i), ilut)

            pos = binary_search(ilut_list, ilut, nifd+1)

            if (pos > 0) then 
                call extract_sign(ilut_list(:,pos), temp_sign)

                nI_weights(i) = temp_sign(1)

            end if
        end do

    end subroutine find_states_in_list

    subroutine get_highest_pop(n_states, largest_dets, norm)
        ! routine to give the n_states most populated states largest_dets
        ! globally
        integer, intent(in) :: n_states
        integer(n_int), intent(out) :: largest_dets(0:niftot, n_states)
        real(dp), intent(out), optional :: norm
#ifdef __DEBUG
        character(*), parameter :: this_routine = "get_highest_pop"
#endif
        integer(n_int) :: largest_dets_node(0:niftot,n_states)
        real(dp) :: norm_node

        call get_highest_pop_node(n_states, largest_dets_node, norm_node)

        call find_highest_sign_per_node(n_states, largest_dets_node, largest_dets)

    end subroutine get_highest_pop

    subroutine find_highest_sign_per_node(n_states, largest_dets_node, largest_dets)
        ! routine to find the largest signs on each node and store them 
        ! sequentially into the global list 
        integer, intent(in) :: n_states
        integer(n_int), intent(inout) :: largest_dets_node(0:niftot,n_states)
        integer(n_int), intent(out) :: largest_dets(0:niftot, n_states)

        integer :: i, max_pos, j
        real(dp) :: tmp_sign(lenof_sign), max_sign
        real(dp) :: reduce_in(2), reduce_out(2)
        integer(n_int) :: max_det(0:niftot)


        do i = 1, n_states
            max_sign = 0.0_dp
            max_pos = 1 

            do j = n_states, 1, -1
                call extract_sign(largest_dets_node(:,j), tmp_sign)
                ! why is this call?
                if (any(largest_dets_node(:,j) /= 0)) then
 
#ifdef __CMPLX
                    max_sign = sqrt(sum(abs(tmp_sign(1::2)))**2 + sum(abs(tmp_sign(2::2)))**2)
#else
                    max_sign = sum(real(abs(tmp_sign),dp))
#endif

                    ! We have the largest sign
                    max_pos = j
                    exit
                end if
            end do
            reduce_in = [max_sign, real(iProcIndex,dp)]
            call MPIAllReduceDatatype(reduce_in,1,MPI_MAXLOC,MPI_2DOUBLE_PRECISION,reduce_out)

            if (iProcIndex == nint(reduce_out(2))) then
                max_det = largest_dets_node(:,max_pos)
                ! and then set it to zero
                largest_dets_node(:,max_pos) = 0
            else 
                max_det = 0
            end if

            call MPIBCast(max_det ,NIfTot+1,nint(reduce_out(2)))

            if (iProcIndex == root) then 
                largest_dets(:,i) = max_det
            end if
        end do

    end Subroutine find_highest_sign_per_node

    subroutine get_highest_pop_node(n_states, largest_dets, all_norm)
        ! routine to give the n_states most populated states largest_dets
        ! per node
        integer, intent(in) :: n_states
        integer(n_int), intent(out) :: largest_dets(0:niftot, n_states)
        real(dp), intent(out), optional :: all_norm
#ifdef __DEBUG
        character(*), parameter :: this_routine = "get_highest_pop_node"
#endif
        real(dp) :: norm

        call return_most_populated_states(n_symmetry_states, &
            largest_dets, CurrentDets, TotWalkers, norm)

        if (present(all_norm)) then
            call MpiSum(norm, all_norm)
        end if

    end subroutine get_highest_pop_node

    function apply_rotation(in_orbs) result(out_orbs)
        ! function to rotate the k- or r-vectors and return the mapped 
        ! orbitals 
        integer, intent(in) :: in_orbs(nBasis/2)
        integer :: out_orbs(nBasis/2)
#ifdef __DEBUG
        character(*), parameter :: this_routine = "apply_rotation"
#endif
        integer :: i

        if (.not. t_symmetry_rotation .or. (symmetry_rotation_angle == 0.0_dp & 
            .or. symmetry_rotation_angle == 360.0_dp)) then 
            out_orbs = in_orbs
            return
        end if

        ! i need the k- or r-vectors (or provide them as input? TDB)
        do i = 1, nBasis/2
            out_orbs(i) = rotate(in_orbs(i))
        end do

    end function apply_rotation

    function rotate_orb(in_orb) result(out_orb)
        ! function to actually apply the rotation to the basis vectors
        integer, intent(in) :: in_orb
        integer :: out_orb
#ifdef __DEBUG
        character(*), parameter :: this_routine = "rotate"
#endif
        integer :: vec(3), rot_vec(3)

        ASSERT(associated(lat))

        if (lat%is_k_space()) then 
            vec = lat%get_k_vec(in_orb)
        else
            vec = lat%get_r_vec(in_orb)
        end if

        rot_vec = rotate(vec)

        ! apply pbc (should be done within the get_orb_from_k_vec function i hope..
        out_orb = lat%get_orb_from_k_vec(rot_vec)

    end function rotate_orb

    function rotate_vec(in_vec) result(out_vec)
        integer, intent(in) :: in_vec(3)
        integer :: out_vec(3)
        real(dp) :: rot_mat(2,2), angle

        angle = symmetry_rotation_angle * pi / 180.0

        rot_mat(1,1) = cos(angle)
        rot_mat(1,2) = -sin(angle)
        rot_mat(2,1) = sin(angle)
        rot_mat(2,2) = cos(angle)

        ! apply the actual rotation to the vector..
        out_vec(1:2) = nint(matmul(rot_mat, real(in_vec(1:2))))
        out_vec(3) = 0

    end function rotate_vec

end module analyse_wf_symmetry

