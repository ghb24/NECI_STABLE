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

    use bit_rep_data, only: niftot, nifd, noffsgn, nifsgn

    use constants, only: n_int, dp, pi, lenof_sign

    use util_mod, only: binary_search, binary_search_int

    use bit_reps, only: extract_sign, encode_sign, decode_bit_det

    use DetBitOps, only: EncodeBitDet, ilut_lt, ilut_gt, DetBitEq

    use sort_mod, only: sort

    use unit_test_helpers, only: print_matrix

    use ras, only: sort_orbitals

    use hist, only: ssquared_contrib

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

    real(dp) :: mirror_x(2,2) = reshape([1.0,0.0,0.0,-1.0],[2,2])
    real(dp) :: mirror_y(2,2) = reshape([-1.0,0.0,0.0,1.0],[2,2])
    real(dp) :: mirror_d(2,2) = reshape([0.0,-1.0,-1.0,0.0],[2,2])
    real(dp) :: mirror_o(2,2) = reshape([0.0,1.0,1.0,0.0],[2,2])

    real(dp) :: inv_matrix(2,2) = reshape([-1.0,0.0,0.0,-1.0],[2,2])

    interface inversion
        module procedure inversion_orb
        module procedure inversion_vec
    end interface inversion

    interface rotate
        module procedure rotate_orb
        module procedure rotate_vec
    end interface rotate

    interface mirror
        module procedure :: mirror_orb
        module procedure :: mirror_vec
    end interface mirror

contains

    subroutine print_point_group_matrix_rep(states)
        ! output the symmetry operation of a specific point group in a given
        ! basis
        integer, intent(in) :: states(:,:)
        character(*), parameter :: this_routine = "print_point_group_matrix_rep"

        if (lat%get_ndim() == 2) then 
            ! for 2D it is the 8 fold point group symmetry and S^2 EV..
            ! activate all the symmetries
            t_symmetry_mirror = .true.
            t_symmetry_rotation = .true.
            t_symmetry_inversion = .true.

            call print_d4h_pg(states)
        else
            call stop_all(this_routine, "not yet implemented!")
        end if


    end subroutine print_point_group_matrix_rep

    subroutine print_d4h_pg(states)
        ! construct the d4h symmetry operation matrix representations
        ! the operation for now are: E, 2C4, C2, Mv, Md, and inversion
        integer, intent(in) :: states(:,:)

        integer :: i, orig_orbs(nBasis/2), trans_orbs(nBasis/2), &
                   matrix_rep(size(states,2),size(states,2)), & 
                   temp_states(nel,size(states,2)), phase

        orig_orbs = get_spatial(brr(1:nBasis:2))

        ! first E:
        trans_orbs = orig_orbs
        matrix_rep = construct_matrix_representation(states, orig_orbs, trans_orbs)
        print *, "E:"
        call print_matrix(matrix_rep)
        print *, "character: ", sum([(matrix_rep(i,i), i = 1, size(matrix_rep,1))])

        ! first C4:
        trans_orbs = apply_rotation(orig_orbs, 90.0_dp)
        matrix_rep = construct_matrix_representation(states, orig_orbs, trans_orbs)
        print *, "C4: "
        call print_matrix(matrix_rep)
        print *, "character: ", sum([(matrix_rep(i,i), i = 1, size(matrix_rep,1))])

!         ! first C4:
!         trans_orbs = apply_rotation(orig_orbs, 270.0)
!         matrix_rep = construct_matrix_representation(states, orig_orbs, trans_orbs)
!         print *, "2C4: "
!         call print_matrix(matrix_rep)
!         print *, "character: ", sum([(matrix_rep(i,i), i = 1, size(matrix_rep,1))])

        ! C2:
        trans_orbs = apply_rotation(orig_orbs, 180.0_dp)
        matrix_rep = construct_matrix_representation(states, orig_orbs, trans_orbs)
        print *, "C2: "
        call print_matrix(matrix_rep)
        print *, "character: ", sum([(matrix_rep(i,i), i = 1, size(matrix_rep,1))])

        ! Mv
        trans_orbs = apply_mirror(orig_orbs, 'x')
        matrix_rep = construct_matrix_representation(states, orig_orbs, trans_orbs)
        print *, "Mx: "
        call print_matrix(matrix_rep)
        print *, "character: ", sum([(matrix_rep(i,i), i = 1, size(matrix_rep,1))])

        ! Md
        trans_orbs = apply_mirror(orig_orbs, 'd')
        matrix_rep = construct_matrix_representation(states, orig_orbs, trans_orbs)
        print *, "Md: "
        call print_matrix(matrix_rep)
        print *, "character: ", sum([(matrix_rep(i,i), i = 1, size(matrix_rep,1))])
        
        ! i
        trans_orbs = apply_inversion(orig_orbs)
        matrix_rep = construct_matrix_representation(states, orig_orbs, trans_orbs)
        print *, "i: "
        call print_matrix(matrix_rep)
        print *, "character: ", sum([(matrix_rep(i,i), i = 1, size(matrix_rep,1))])

    end subroutine print_d4h_pg

    function construct_matrix_representation(states, orig_orbs, trans_orbs) & 
            result(matrix)
        ! construct the matrix representation of the symmetry operation 
        ! for a given basis. the symmetry operation is encoded in 
        ! orig orbs and trans orbs
        integer, intent(in) :: states(:,:), orig_orbs(nBasis/2), trans_orbs(nBasis/2)
        integer :: matrix(size(states,2),size(states,2))
        integer :: i, j, nJ(nel), phase

        matrix = 0

        do i = 1, size(states,2)
            do j = 1, size(states,2)
                call apply_transformation(states(:,j), orig_orbs, trans_orbs, &
                    nJ, phase)

                if (all(nJ == states(:,i))) then
                    matrix(i,j) = phase
!                     matrix(j,i) = phase
                end if
            end do
        end do
        
    end function construct_matrix_representation

    subroutine analyze_full_wavefunction_sym(sym_labels, ilut_list_opt)
        ! give the symmetry eigenvalues of a wavefunction in ilut_format 
        ! of certain symmetry operations, depending on the lattice point group
        ! e.g. for 2D: rot90, rot180, rot270, m_x, m_y, m_d, m_o and the S^2 EV
        ! on can either provide a wavfunction in ilut_list, or otherwise 
        ! it is obtained from the FCIQMC wavefunction depending on the 
        ! input
        real(dp), intent(out), allocatable :: sym_labels(:)
        integer(n_int), intent(in), optional :: ilut_list_opt(:,:)
        character(*), parameter :: this_routine = "analyze_full_wavefunction_sym"
        integer(n_int), allocatable :: ilut_list(:,:), trans_pg_wf(:,:,:), &
                                       ilut_spin(:,:)
        integer :: n_syms, i, ms, nI(nel)

        ASSERT(associated(lat))

        if (lat%get_ndim() == 2) then 
            ! for 2D it is the 8 fold point group symmetry and S^2 EV..
            n_syms = 10
            ! activate all the symmetries
            t_symmetry_mirror = .true.
            t_symmetry_rotation = .true.
            t_symmetry_inversion = .true.
        else
            call stop_all(this_routine, "not yet implemented!")
        end if

        allocate(sym_labels(n_syms))
        sym_labels = 0.0_dp

        if (present(ilut_list_opt)) then 
            allocate(ilut_list(0:niftot, size(ilut_list_opt,2)), &
                source = ilut_list_opt)

        else
            call init_symmetry_states()
            allocate(ilut_list(0:NIfTot, size(symmetry_states_ilut,2)), &
                source = symmetry_states_ilut)

        end if

        if (lat%get_ndim() == 2) then 
            trans_pg_wf = apply_2D_point_group(ilut_list)
        else
            call stop_all(this_routine, "not yet implemented!")
        end if

        ! now we have to calculate <y|y'> to get the symmetry EV
        ! but dont do the S^2 yet..
        do i = 1, n_syms - 1
            sym_labels(i) = calc_overlap(ilut_list, trans_pg_wf(:,:,i))
        end do

        do i = 1, size(ilut_list,2)
            sym_labels(n_syms) = sym_labels(n_syms) + &
                ssquared_contrib(ilut_list(:,i), .false., size(ilut_list,2), &
                ilut_list)
        end do

        ! and we need the S_z^2 contribution: 
        call decode_bit_det(nI, ilut_list(:,1))
        ms = sum(get_spin_pn(nI))
        sym_labels(n_syms) = sym_labels(n_syms) + real(ms * (ms + 2) / 4.0, dp)
            
!         ilut_spin = apply_s_squared(ilut_list)
! 
!         sym_labels(n_syms) = calc_overlap(ilut_list, ilut_spin)

    end subroutine analyze_full_wavefunction_sym

    function apply_s_squared(ilut_list) result(ilut_spin)
        ! function to apply the S^2 operator to a given wavefunction 
        ! figure that out, what we have to do here.. 
        integer(n_int), intent(inout) :: ilut_list(:,:)
        integer(n_int), allocatable :: ilut_spin(:,:)

        allocate(ilut_spin(0:niftot, size(ilut_list,2)))
        ilut_spin = 0_n_int

    end function apply_s_squared

    function calc_overlap(ilutI, ilutJ) result(overlap)
        ! calculate the overlap between two wavefunction I and J
        integer(n_int), intent(in) :: ilutI(:,:), ilutJ(:,:)
        real(dp) :: overlap 
        real(dp) :: signI(lenof_sign), signJ(lenof_sign)

        integer :: i, pos
        ! i am pretty sure the lists are ordered so I can binary search.. 

        overlap = 0.0_dp

        do i = 1, size(ilutI,2)
            pos = binary_search(ilutJ, ilutI(:,i), nifd+1)

            if (pos > 0) then 
                call extract_sign(ilutI(:,i), signI)
                call extract_sign(ilutJ(:,pos), signJ)

                overlap = overlap + signI(1)*signJ(1)

            end if
        end do

    end function calc_overlap

    function apply_2D_point_group(ilut_list) result(trans_wf)
        ! apply the point group symmetries of a 2D square lattice
        integer(n_int), intent(inout) :: ilut_list(:,:)
        integer(n_int), allocatable :: trans_wf(:,:,:)

        integer :: i, j
        real(dp) :: rot_angle
        character(4) :: mir_axis
        integer :: sort_ind(size(ilut_list,2))
        integer :: matrix(size(ilut_list,2),size(ilut_list,2))
        real(dp) :: signI(lenof_sign), signJ(lenof_sign)

        allocate(trans_wf(0:niftot,size(ilut_list,2), 9))

        trans_wf = 0_n_int

        rot_angle = 0.0_dp

        do i = 1, 4 
            trans_wf(:,:,i) = apply_rotation_wf(ilut_list, rot_angle)

            rot_angle = rot_angle + 90.0_dp

        end do

        mir_axis = 'xydo'

        do i = 5, 8
            trans_wf(:,:,i) = apply_mirror_wf(ilut_list, mir_axis(i-4:i-4), sort_ind)
        end do

        trans_wf(:,:,9) = apply_inversion_wf(ilut_list)

        ! test to output the matrix representation of the symmetry operation
!         print *, "sort_ind: ", sort_ind
!         matrix = 0
!         do i = 1, size(ilut_list,2)
!             call extract_sign(ilut_list(:,i), signI)
!             call extract_sign(trans_wf(:,sort_ind(i),8), signJ)
!             
! !             matrix(i,sort_ind(i)) = int(sign(1.0_dp, signI(1)*signJ(1)))
! 
!         end do


!         print *, "m_o: "
!         call print_matrix(matrix)


    end function apply_2D_point_group

    function apply_inversion_wf(ilut_list, sort_ind) result(inv_wf)
        ! apply inversion symmetry to a wavefunction 
        integer(n_int), intent(inout) :: ilut_list(:,:)
        integer, intent(out), optional :: sort_ind(size(ilut_list,2))
        integer(n_int) :: inv_wf(0:size(ilut_list,1)-1,size(ilut_list,2))

        integer :: orig_orbs(nBasis/2), trans_orbs(nBasis/2), i, &
                   orig_states(nel,size(ilut_list,2)), inv_states(nel, size(ilut_list,2))
        real(dp) :: orig_weights(size(ilut_list,2)), temp_sign(lenof_sign), &
                    inv_weights(size(ilut_list,2))


        orig_orbs = get_spatial(brr(1:nBasis:2))

        trans_orbs = apply_inversion(orig_orbs)

        ! decode the original information
        do i = 1, size(ilut_list,2)
            call decode_bit_det(orig_states(:,i), ilut_list(:,i))
            call extract_sign(ilut_list(:,i), temp_sign)
            orig_weights(i) = temp_sign(1)
        end do

        if (present(sort_ind)) then
            call transform_states(orig_orbs, trans_orbs, size(ilut_list,2), &
                orig_states, orig_weights, ilut_list, inv_states, inv_weights, &
                inv_wf, sort_ind)
        else
            call transform_states(orig_orbs, trans_orbs, size(ilut_list,2), &
                orig_states, orig_weights, ilut_list, inv_states, inv_weights, inv_wf)
        end if


    end function apply_inversion_wf

    function apply_mirror_wf(ilut_list, mirror_axis, sort_ind) result(mir_wf)
        ! function to apply a mirror symmetry to the given wavefunction 
        ! encoded in ilut_list
        integer(n_int), intent(inout) :: ilut_list(:,:)
        character(1), intent(in) :: mirror_axis
        integer, intent(out), optional :: sort_ind(size(ilut_list,2))
        integer(n_int) :: mir_wf(size(ilut_list,1),size(ilut_list,2))
#ifdef __DEBUG
        character(*), parameter :: this_routine = "apply_mirror_wf"
#endif
        integer :: orig_orbs(nBasis/2), trans_orbs(nBasis/2), i, &
                   orig_states(nel,size(ilut_list,2)), mir_states(nel, size(ilut_list,2))
        real(dp) :: orig_weights(size(ilut_list,2)), temp_sign(lenof_sign), &
                    mir_weights(size(ilut_list,2))

        orig_orbs = get_spatial(brr(1:nBasis:2))

        trans_orbs = apply_mirror(orig_orbs, mirror_axis)

        ! decode the original information
        do i = 1, size(ilut_list,2)
            call decode_bit_det(orig_states(:,i), ilut_list(:,i))
            call extract_sign(ilut_list(:,i), temp_sign)
            orig_weights(i) = temp_sign(1)
        end do

        if (present(sort_ind)) then
            call transform_states(orig_orbs, trans_orbs, size(ilut_list,2), &
                orig_states, orig_weights, ilut_list, mir_states, mir_weights, &
                mir_wf, sort_ind)
        else
            call transform_states(orig_orbs, trans_orbs, size(ilut_list,2), &
                orig_states, orig_weights, ilut_list, mir_states, mir_weights, mir_wf)
        end if

    end function apply_mirror_wf

    function apply_rotation_wf(ilut_list, rot_angle) result(rot_wf)
        ! function to apply a rotation of given angle to the whole wavefunction
        ! encoded in ilut_list
        integer(n_int), intent(inout) :: ilut_list(:,:)
        real(dp), intent(in) :: rot_angle
        integer(n_int) :: rot_wf(size(ilut_list,1),size(ilut_list,2))
#ifdef __DEBUG
        character(*), parameter :: this_routine = "apply_rotation_wf"
#endif
        integer :: orig_orbs(nBasis/2), trans_orbs(nBasis/2), i, &
                   orig_states(nel,size(ilut_list,2)), rot_states(nel, size(ilut_list,2))
        real(dp) :: orig_weights(size(ilut_list,2)), temp_sign(lenof_sign), &
                    rot_weights(size(ilut_list,2))

        orig_orbs = get_spatial(brr(1:nBasis:2))

        trans_orbs = apply_rotation(orig_orbs, rot_angle)

        ! decode the original information
        do i = 1, size(ilut_list,2)
            call decode_bit_det(orig_states(:,i), ilut_list(:,i))
            call extract_sign(ilut_list(:,i), temp_sign)
            orig_weights(i) = temp_sign(1)
        end do

        call transform_states(orig_orbs, trans_orbs, size(ilut_list,2), &
            orig_states, orig_weights, ilut_list, rot_states, rot_weights, rot_wf)

    end function apply_rotation_wf

    function rot_matrix(rot_angle) result(mat)
        real(dp), intent(in) :: rot_angle
        real(dp) :: mat(2,2)

        real(dp) :: angle

        angle = rot_angle * pi / 180.0
        mat = reshape([cos(angle),sin(angle),-sin(angle),cos(angle)],[2,2])
!         mat(1,1) = cos(angle)
!         mat(1,2) = -sin(angle)
!         mat(2,1) = sin(angle)
!         mat(2,2) = cos(angle)

    end function rot_matrix
        
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
            sym_orbs = apply_rotation(orig_orbs, symmetry_rotation_angle)
            
            print *, "orig orbs -> sym_orbs: "
            do i = 1, nBasis/2
                print *, orig_orbs(i), " -> ", sym_orbs(i)
            end do

            ! i can ofc apply multiple symmetries
            ! MIRROR
            sym_orbs = apply_mirror(sym_orbs, symmertry_mirror_axis)

            ! INVERSION
            sym_orbs = apply_inversion(sym_orbs)

            ! now i have the mapping between the original and the 
            ! symmetry transformed orbitals

            call transform_states(orig_orbs, sym_orbs, n_symmetry_states, & 
                symmetry_states, symmetry_weights, symmetry_states_ilut, transformed_states, &
                transformed_weights, transformed_states_ilut)

            call compare_states(n_symmetry_states, symmetry_states_ilut, transformed_states_ilut)

        end if

    end subroutine analyze_wavefunction_symmetry

    subroutine compare_states(n_states, orig_states, trans_states)
        ! compare the original and the symmetry transformed states 
        integer, intent(in) :: n_states
        integer(n_int), intent(in) :: orig_states(0:niftot, n_states), &
                                      trans_states(0:niftot, n_states)

        integer(n_int) :: null_int(0:niftot), ilutI(0:niftot), ilutJ(0:niftot)
        integer :: i, j, k, l
        ! print the original and the transformed states next to each other 
        ! and missing determinants in the repective lists are indicated 
        i = 1 
        j = 1

        null_int = 0_n_int

        ! or first create a list and an integer indicator, where the 
        ! determinant is from.. 
        do k = 1, 2*n_states
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
            if (i == j .and. i > n_states) exit

            if (i > n_states .and. j <= n_states) then 
                ! if i is already above the list then the rest of the 
                ! entries are from list J
                do l = j, n_states
                    ilutJ = trans_states(:,l)
                    call print_2_states(null_int, ilutJ)
                end do
                exit
            end if
            ! and if j is alreay above, then the rest is from I 
            if (i <= n_states .and. j > n_states) then 
                do l = i, n_states
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

    subroutine transform_states(orig_orbs, transformed_orbs, n_states, orig_states, &
            orig_weights, orig_iluts, transformed_states, transformed_weights, &
            transformed_states_ilut, sort_ind)
        integer, intent(in) :: orig_orbs(nBasis/2), transformed_orbs(nBasis/2)
        integer, intent(in) :: n_states
        integer, intent(inout) :: orig_states(nel,n_states)
        real(dp), intent(inout) :: orig_weights(n_states)
        integer(n_int), intent(inout) :: orig_iluts(0:niftot,n_states)
        integer, intent(out) :: transformed_states(nel, n_states)
        real(dp), intent(out) :: transformed_weights(n_states)
        integer, intent(out), optional :: sort_ind(n_states)
        integer(n_int), intent(out) :: transformed_states_ilut(0:niftot,n_states)
#ifdef __DEBUG
        character(*), parameter :: this_routine = "transform_states"
#endif
        integer :: i, n_phase, ind(n_states)
        real(dp) :: tmp_sign(lenof_sign)
    
        tmp_sign = 0.0_dp

        ! do it plain ans stupid for now.. 
        do i = 1, n_states

             call apply_transformation(orig_states(:,i), orig_orbs, &
                 transformed_orbs, transformed_states(:,i), n_phase)

             transformed_weights(i) = real(n_phase,dp) * orig_weights(i)

             tmp_sign(1) = transformed_weights(i)

             call EncodeBitDet(transformed_states(:,i), transformed_states_ilut(:,i))
             call encode_sign(transformed_states_ilut(:,i), tmp_sign)

        end do

        ! the original highest pop list is sorted by weight i guess.. 
        ! sort them by the integers in the ilut representation
        ind = [(i, i = 1, n_states)]

        call sort(orig_iluts, ind)

        orig_states = orig_states(:,ind)
        orig_weights = orig_weights(ind)


        ind = [(i, i = 1, n_states)]

        call sort(transformed_states_ilut, ind)

        transformed_states = transformed_states(:,ind)
        transformed_weights = transformed_weights(ind)

        if (present(sort_ind)) then
            sort_ind = ind
        end if

        
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
            if (pos <= 0) then 
                print *, "orig_orbs: ", orig_orbs
                print *, "nI: ", nI
                print *, "spatial: ", get_spatial(nI)
            end if
            ASSERT(pos > 0)
            
            if (is_beta(nI(i))) then 
                nJ(i) = 2*trans_orbs(pos) - 1
            else 
                nJ(i) = 2*trans_orbs(pos)
            end if

        end do

        if (present(n_phase)) then 
            call sort_orbitals(nJ, n_phase)
!             call sort(nJ, par = n_phase)
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

    function apply_inversion(in_orbs) result(out_orbs)
        ! apply inversion through the k-point
        integer, intent(in) :: in_orbs(nBasis/2)
        integer :: out_orbs(nBasis/2)

        integer :: i 

        if (.not. t_symmetry_inversion) then 
            out_orbs = in_orbs
            return
        end if

        do i = 1, nBasis/2
            out_orbs(i) = inversion(in_orbs(i))
        end do

    end function apply_inversion

    function apply_mirror(in_orbs, mirror_axis) result(out_orbs)
        ! function to mirror the k- or r- vectors along a specified axis
        integer, intent(in) :: in_orbs(nBasis/2)
        character(1), intent(in) :: mirror_axis
        integer :: out_orbs(nBasis/2)
#ifdef __DEBUG
        character(*), parameter :: this_routine = "apply_mirror"
#endif
        integer :: i

        if (.not. t_symmetry_mirror .or. mirror_axis == '0') then
            out_orbs = in_orbs
            return
        end if

        do i = 1, nBasis/2
            out_orbs(i) = mirror(in_orbs(i), mirror_axis)
        end do

    end function apply_mirror

    function apply_rotation(in_orbs, rot_angle) result(out_orbs)
        ! function to rotate the k- or r-vectors and return the mapped 
        ! orbitals 
        integer, intent(in) :: in_orbs(nBasis/2)
        real(dp), intent(in) :: rot_angle
        integer :: out_orbs(nBasis/2)
#ifdef __DEBUG
        character(*), parameter :: this_routine = "apply_rotation"
#endif
        integer :: i

        if (.not. t_symmetry_rotation .or. (rot_angle == 0.0_dp & 
            .or. rot_angle == 360.0_dp)) then 
            out_orbs = in_orbs
            return
        end if

        ! i need the k- or r-vectors (or provide them as input? TDB)
        do i = 1, nBasis/2
            out_orbs(i) = rotate(in_orbs(i), rot_angle)
        end do

    end function apply_rotation


    function mirror_orb(in_orb, mirror_axis) result(out_orb)
        ! function to actually apply the mirroring to an orbital
        integer, intent(in) :: in_orb
        character(1), intent(in) :: mirror_axis
        integer :: out_orb
#ifdef __DEBUG
        character(*), parameter :: this_routine = "mirror_orb"
#endif
        integer :: vec(3), mir_vec(3)

        ASSERT(associated(lat))

        if (mirror_axis == '0') then
            out_orb = in_orb
            return
        end if

        if (lat%is_k_space()) then
            vec = lat%get_k_vec(in_orb)
        else
            vec = lat%get_r_vec(in_orb)
        end if

        mir_vec = mirror(vec, mirror_axis)

        out_orb = lat%get_orb_from_k_vec(mir_vec)

    end function mirror_orb

    function mirror_vec(in_vec, mirror_axis) result(out_vec)
        integer, intent(in) :: in_vec(3)
        character(1), intent(in) :: mirror_axis
        integer :: out_vec(3)

        select case (mirror_axis)
        case('x')
            out_vec(1:2) = nint(matmul(mirror_x, real(in_vec(1:2))))

        case('y')
            out_vec(1:2) = nint(matmul(mirror_y, real(in_vec(1:2))))

        case('d')
            out_vec(1:2) = nint(matmul(mirror_d, real(in_vec(1:2))))

        case('o')
            out_vec(1:2) = nint(matmul(mirror_o, real(in_vec(1:2))))

        case ('0') 
            out_vec = in_vec

        case Default

            call stop_all("mirror_vec", "incorrect mirroring axis!")

        end select

        out_vec(3) = in_vec(3)

        
    end function mirror_vec

    function inversion_orb(in_orb) result(out_orb)
        integer, intent(in) :: in_orb
        integer :: out_orb
#ifdef __DEBUG
        character(*), parameter :: this_routine = "inversion_orb"
#endif
        integer :: vec(3), inv_vec(3)

        ASSERT(associated(lat))

        if (lat%is_k_space()) then
            vec = lat%get_k_vec(in_orb)
        else
            vec = lat%get_r_vec(in_orb)
        end if

        inv_vec = inversion(vec)

        out_orb = lat%get_orb_from_k_vec(inv_vec)

    end function inversion_orb

    function inversion_vec(in_vec) result(out_vec)
        integer, intent(in) :: in_vec(3)
        integer :: out_vec(3)

        out_vec(1:2) = nint(matmul(inv_matrix, real(in_vec(1:2))))
        out_vec(3) = in_vec(3)

    end function inversion_vec
                
    function rotate_orb(in_orb, rot_angle) result(out_orb)
        ! function to actually apply the rotation to the basis vectors
        integer, intent(in) :: in_orb
        real(dp), intent(in) :: rot_angle
        integer :: out_orb
#ifdef __DEBUG
        character(*), parameter :: this_routine = "rotate_orb"
#endif
        integer :: vec(3), rot_vec(3)

        ASSERT(associated(lat))

        if (rot_angle == 0.0_dp .or. rot_angle == 360.0_dp) then 
            out_orb = in_orb
            return
        end if

        if (lat%is_k_space()) then 
            vec = lat%get_k_vec(in_orb)
        else
            vec = lat%get_r_vec(in_orb)
        end if

        rot_vec = rotate(vec, rot_angle)

        ! apply pbc (should be done within the get_orb_from_k_vec function i hope..
        out_orb = lat%get_orb_from_k_vec(rot_vec)

    end function rotate_orb

    function rotate_vec(in_vec, rot_angle) result(out_vec)
        integer, intent(in) :: in_vec(3)
        real(dp), intent(in) :: rot_angle
        integer :: out_vec(3)

        ! apply the actual rotation to the vector..
        out_vec(1:2) = nint(matmul(rot_matrix(rot_angle), real(in_vec(1:2))))
        out_vec(3) = in_vec(3)

    end function rotate_vec

end module analyse_wf_symmetry

