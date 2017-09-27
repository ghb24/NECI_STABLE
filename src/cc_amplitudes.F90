#include "macros.h" 

module cc_amplitudes 

    use SystemData, only: nbasis, nel, nOccAlpha, nOccBeta, ElecPairs
    use detbitops, only: get_bit_excitmat, FindBitExcitLevel
    use FciMCData, only: totwalkers, ilutref, currentdets, AllNoatHf, projedet, &
                         HashIndex, CurrentDets, ll_node
    use bit_reps, only: extract_sign
    use constants, only: dp, lenof_sign, EPS, n_int, bits_n_int
    use bit_reps, only: niftot, nifdbo, decode_bit_det
    use replica_data, only: AllEXLEVEL_WNorm
    use back_spawn, only: setup_virtual_mask, mask_virt_ni
    use hash, only: hash_table_lookup, FindWalkerHash
    use util_mod, only: swap

    implicit none 

    logical :: t_cc_amplitudes = .false. 
    logical :: t_plot_cc_amplitudes = .false.

    ! we need the order of cluster operators..
    integer :: cc_order = 0 

    integer :: cc_delay = 1000

    integer :: est_triples = 0, est_quads = 0, est_doubles = 0
    integer :: n_singles = 0, n_doubles = 0
    integer :: n_triples = 0, n_quads = 0

    logical :: t_store_hash_quadrupels = .true. 
    integer :: quad_hash_size

    ! i guess i need a new hash table for cc amplitudes, since i also 
    ! want the amplitudes.. but how do i deal with clashes here? 
    type cc_hash
        logical :: found = .false.
        real(dp) :: amp = 0.0_dp
        integer(n_int), allocatable :: ind(:)
        type(cc_hash), pointer :: next => null()
    end type cc_hash

    type(cc_hash), pointer :: quad_hash(:)

    integer :: n_clashes = 0

!     interface get_amp 
!         module procedure get_amp_ind
!         module procedure get_amp_ex
!     end interface
 
    ! maybe it would be nice to have a type which encodes this information 
    ! and which gives an easy and nice way to de/encode the indices 
    ! involved.. 
    type cc_amplitude 
        integer :: order = 0
        integer, allocatable :: operators(:,:,:)
        real(dp), allocatable :: amplitudes(:)
        logical, allocatable :: set_flag(:)

        integer :: n_ops = 0

    contains 

        procedure :: get_ex
        procedure :: get_ind
        procedure :: get_amp => get_amp_ind

    end type cc_amplitude

    ! and make a global cc_ops 
    type(cc_amplitude), allocatable :: cc_ops(:)

    integer, allocatable :: ind_matrix_singles(:,:)
    integer, allocatable :: ind_matrix_doubles(:,:) 

    integer :: n_virt_pairs = 0, nvirt = 0

    integer, allocatable :: elec_ind_mat(:,:), orb_ind_mat(:,:)

contains 

    subroutine print_cc_amplitudes 
        ! routine to print out the cc-amplitudes to check there 
        ! magnitude 
        use util_mod, only: get_free_unit, get_unique_filename
        integer :: i, j, iunit
        character(12) :: filename
        character(1) :: x1
        type(cc_hash), pointer :: temp_node
        
        if (.not. allocated(cc_ops)) then 
            print *, "cc amplitudes not yet allocated! cant print!"
            return
        end if

        do i = 1, 3
            ! open 4 files to print the cc-amps
            iunit = get_free_unit() 
            write(x1,'(I1)') i 
            call get_unique_filename('cc_amps_'//trim(x1), .true., .true., 1, &
                filename)
            open(iunit, file = filename, status = 'unknown')

            do j = 1, cc_ops(i)%n_ops
                if (abs(cc_ops(i)%get_amp(j)) < EPS) cycle
                write(iunit, '(f15.7)') abs(cc_ops(i)%get_amp(j))
            end do

            close(iunit)

        end do

        iunit = get_free_unit() 
        call get_unique_filename('cc_amps_4', .true., .true., 1, &
            filename)

        open(iunit, file = filename, status = 'unknown')

        ! i have to do the quads special since it is stored in a hash.. 
        print *, "hash size: ", quad_hash_size
        do i = 1, quad_hash_size
            temp_node => quad_hash(i) 
            
            if (temp_node%found) then 
                write(iunit, '(f15.7)') abs(temp_node%amp)
            
            end if 

            do while(associated(temp_node%next)) 
                if (temp_node%next%found) then 
                    write(iunit, '(f15.7)') abs(temp_node%next%amp)
                end if

                temp_node => temp_node%next

            end do
        end do
        close(iunit)
        
    end subroutine print_cc_amplitudes

    function cc_singles_factor() result(factor) 
        ! this function should provide the correct factor to the 
        ! cepa-shift for the singles.. if the correct variable are not 
        ! yet set or sampled as 0 it should default to 1
        real(dp) :: factor
        character(*), parameter :: this_routine = "cc_singles_factor" 

        real(dp) :: fac_triples, fac_doubles, weight 

        weight = 1.0_dp / 4.0_dp

        ! the singles should be influenced by the triples and doubles.. 
        ! but this i have not figured out correctly.. 
        ! so for now return 1 always 
        factor = 1.0_dp 
        return

        if (est_triples == 0) then 
            fac_triples = 0.0_dp 

            ! fix the weight in this case
            weight = 1.0_dp

        else 
            ! for now just deal with the L^0 norm of the triples 
            fac_triples = min(AllEXLEVEL_WNorm(0,3,1) / real(est_triples, dp), 1.0_dp)

        end if

        ! with the doubles i am scared that the estimated number of 
        ! doubles could actually be lower then the sampled ones.. 
        if (est_doubles == 0) then

            fac_doubles = 0.0_dp

            weight = 0.0_dp

        else
            ! but 1 should be the maximum..
            fac_doubles = min(AllEXLEVEL_WNorm(0,2,1) / real(est_doubles, dp), 1.0_dp)
        end if

        ! and then we have to combine the two factors with some weighting
        factor = 1.0_dp - (weight * fac_doubles + (1.0_dp - weight)*fac_triples)

    end function cc_singles_factor

    function cc_doubles_factor() result(factor) 
        real(dp) :: factor
        character(*), parameter :: this_routine = "cc_doubles_factor"

        real(dp) :: fac_triples, fac_quads, weight

        weight = 1.0_dp/4.0_dp 

        ! essentiall we would need the triples influence too.. 
        ! but Manu said this is not necessary.. hm.. lets try 
        if (est_triples == 0) then 
            fac_triples = 0.0_dp 
            weight = 0.0_dp 
        else
            fac_triples = min(AllEXLEVEL_WNorm(0,3,1) / real(est_triples, dp), 1.0_dp)
        end if

        if (est_quads == 0) then 
            fac_quads = 0.0_dp 
            weight = 1.0_dp 
        else
            fac_quads = min(AllEXLEVEL_WNorm(0,4,1) / real(est_quads, dp), 1.0_dp)

        end if

        factor = 1.0_dp - (weight * fac_triples + (1.0_dp - weight)*fac_quads)

    end function cc_doubles_factor

    function cc_triples_factor() result(factor)
        real(dp) :: factor 
        character(*), parameter :: this_routine = "cc_triples_factor"

        ! for now just set that to 1 
        factor = 1.0_dp
        
        ! essentially we would need to have the influence of the 
        ! quadrupels and the quintuples.. but those numbers are hard to 
        ! estimate.. 
        ! todo

    end function cc_triples_factor

    function cc_quads_factor() result(factor) 
        real(dp) :: factor
        character(*), parameter :: this_routine = "cc_quads_factor" 

        ! see above
        factor = 1.0_dp

    end function cc_quads_factor 

    subroutine setup_ind_matrix_doubles()
        character(*), parameter :: this_routine = "setup_ind_matrix_doubles"

        integer :: n_elec_pairs, i, j, a, b, elec_i, elec_j
        integer :: ij, ab, orb_a, orb_b, k
        logical :: t_par
        integer(n_int) :: temp_ilut(0:niftot)

        ASSERT(allocated(projedet))
        ASSERT(allocated(iLutRef))
        
        nvirt = nbasis - nel 
        n_virt_pairs = nvirt*(nvirt - 1) / 2
        n_elec_pairs = nel*(nel - 1) / 2 

        call setup_virtual_mask() 

        call setup_elec_ind_mat()
        call setup_orb_ind_mat()

        ! encode only the possible excitations from the closed shell 
        ! reference! 
        ! i < j; a < b and (i,j) /= (a,b) 

        ASSERT(n_elec_pairs == ElecPairs)

        allocate(ind_matrix_doubles(n_elec_pairs,n_virt_pairs))
        ind_matrix_doubles = 0

        temp_ilut = iLutRef(:,1)
        k =  1

        ! i could also just loop over the reference determinant.. 
        do i = 1, nel 
            elec_i = projedet(i,1)
            do j = i + 1, nel 
                elec_j = projedet(j,1)

                t_par = same_spin(elec_i, elec_j)

                ! i have to convert (i,j) -> to a linear index
                ! and i definetly have to use the actual spin orbitals, 
                ! since this is the quantity we have access to in the rest of 
                ! the calculation
                ij = linear_elec_ind(elec_i,elec_j) 
                ! maybe here i could to a check if ij = 0 ?

                ASSERT(ij > 0) 
                ASSERT(ij <= ElecPairs)

                ! now i have to check if the orbitals fit and if the spin 
                ! fits! 
                ! use the virtual mask from back-spawn! 
                if (t_par) then
                    do a = 1, Nvirt
                        orb_a = mask_virt_ni(a,1)
                        do b = a + 1, Nvirt
                            orb_b = mask_virt_ni(b,1)

                            if (same_spin(orb_a, orb_b).and. same_spin(orb_a,elec_i)) then 
                                ab = linear_orb_ind(orb_a,orb_b) 
                                ASSERT(ab > 0)
                                ASSERT(ab <= n_virt_pairs)

                                ind_matrix_doubles(ij,ab) = k 
                                k = k + 1

                            end if
                        end do
                    end do
                else
                    do a = 1, nvirt 
                        orb_a = mask_virt_ni(a,1)
                        do b = a + 1, nvirt 
                            orb_b = mask_virt_ni(b,1) 
                            if (.not. same_spin(orb_a, orb_b)) then 
                                ab = linear_orb_ind(orb_a, orb_b)

                                ASSERT(ab > 0)
                                ASSERT(ab <= n_virt_pairs)

                                ind_matrix_doubles(ij,ab) = k 
                                k = k + 1
                            end if
                        end do
                    end do
                end if
            end do
        end do

        print *, "max double index: ", k - 1

    end subroutine setup_ind_matrix_doubles

    function linear_elec_ind(i,j) result(ind)
        integer, intent(in) :: i, j
        integer :: ind 
        character(*), parameter :: this_routine = "linear_elec_ind" 

        ! for a general orbital indexing i want to get a linear index 
        ! for the electron pairs in the closed shell reference 
        ! i could set up a matrix (nbasis, nbasis) for these pairs 
        ! or i could make a search in the reference det at which position 
        ! (i) and (j) are and use this then to encode the infomation.. 
        ! because i am lazy i think i will setup the matrix! 
        ind = elec_ind_mat(i,j) 

    end function linear_elec_ind

    function linear_orb_ind(a,b) result(ind) 
        integer, intent(in) :: a, b 
        integer :: ind 
        character(*), parameter :: this_routine = "linear_orb_ind" 

        ! see above! 
        ind = orb_ind_mat(a,b) 

    end function linear_orb_ind

    subroutine setup_elec_ind_mat
        character(*), parameter :: this_routine = "setup_elec_ind_mat" 

        integer :: i, j, k, elec_i, elec_j
        if (allocated(elec_ind_mat)) deallocate(elec_ind_mat) 
        ASSERT(allocated(projedet))

        allocate(elec_ind_mat(nbasis,nbasis)) 
        elec_ind_mat = 0

        k = 1 

        do i = 1, nel 
            elec_i = projedet(i,1)
            do j = i + 1, nel 
                elec_j = projedet(j,1)

                elec_ind_mat(elec_i,elec_j) = k 

                k = k + 1
            end do
        end do

        ASSERT(k - 1 == ElecPairs)

    end subroutine setup_elec_ind_mat

    subroutine setup_orb_ind_mat
        character(*), parameter :: this_routine = "setup_orb_ind_mat"

        integer :: i, j, k, orb_i, orb_j 

        if (allocated(orb_ind_mat)) deallocate(orb_ind_mat)

        ASSERT(allocated(mask_virt_ni))

        allocate(orb_ind_mat(nbasis,nbasis)) 
        orb_ind_mat = 0

        k = 1 
        do i = 1, Nvirt
            orb_i = mask_virt_ni(i,1)
            do j = i + 1, nvirt
                orb_j = mask_virt_ni(j,1)

                orb_ind_mat(orb_i, orb_j) = k 
                k = k + 1 
            end do
        end do

        ASSERT(k - 1 == n_virt_pairs)

    end subroutine setup_orb_ind_mat

    function get_amp_ind(this, ind) result(amp)
        ! write an amplitude getter, which gives 0 if the index does not 
        ! fit
        class(cc_amplitude) :: this 
        integer, intent(in) :: ind 
        real(dp) :: amp 

        if (ind == 0) then 
            amp = 0.0_dp
        else 
            amp = this%amplitudes(ind)
        end if

    end function get_amp_ind

    function get_amp_ex(this, elec_ind, orb_ind) result(ind)
        class(cc_amplitude) :: this
        integer, intent(in) :: elec_ind(this%order), orb_ind(this%order)
        real(dp) :: amp

        integer :: ind

        ind = this%get_ind(elec_ind, orb_ind)

        if (ind == 0) then 
            amp = 0.0_dp 
        else
            amp = this%amplitudes(ind)
        end if

    end function get_amp_ex

    subroutine setup_ind_matrix_singles() 
        character(*), parameter :: this_routine = "setup_ind_matrix_singles"

        integer :: i, j, k
        integer(n_int) :: temp_ilut(0:niftot)

        ASSERT(allocated(projedet))
        ASSERT(allocated(iLutRef))

        allocate(ind_matrix_singles(nbasis, nbasis))
        ind_matrix_singles = 0

        temp_ilut = iLutRef(:,1)
        k =  1
        do i = 1, nbasis
            if (IsOcc(temp_ilut, i)) then 
                ! the i is an electron in the reference! 
                ! and for each electron index the virtuals 
                do j = 1, nbasis
                    ! and i should and can specify the spin here!
                    if (IsNotOcc(temp_ilut,j) .and. same_spin(i,j)) then 
                        ind_matrix_singles(i,j) = k 
                        k = k + 1
                    end if
                end do
            end if
        end do

        print *, "max single index: ", k - 1

    end subroutine setup_ind_matrix_singles

    ! and now go to the routines to calculate the number of triples and 
    ! quadrupels 
    subroutine init_cc_amplitudes

        integer :: n_excits(4)

        print *, "testing the cc-amplitudes: "
        print *, "------ test on the cc ------- "

        n_excits = calc_number_of_excitations(nOccAlpha, nOccBeta, 4, & 
            nbasis/2)

        print *, "total number of possible excitations: " 
        print *, n_excits

        call setup_ind_matrix_singles()
        call setup_ind_matrix_doubles()


        call fill_cc_amplitudes()
        call calc_n_triples()
        if (t_store_hash_quadrupels) then 
            quad_hash_size = n_doubles*(n_doubles - 1) / 2
            call init_cc_hash(quad_hash, quad_hash_size)
        end if
        call calc_n_quads()

        print *, "number of hash clashes: ", n_clashes

        print *, "sampled singles: ", n_singles
        print *, "sampled doubles: ", n_doubles

        print *, "est. triples: ", est_triples
        print *, "est. quads:   ", est_quads
        print *, "act. singles: ", AllEXLEVEL_WNorm(0,1,1)
        print *, "act. doubles: ", AllEXLEVEL_WNorm(0,2,1)
        print *, "act. triples: ", AllEXLEVEL_WNorm(0,3,1)
        print *, "act. quads:   ", AllEXLEVEL_WNorm(0,4,1)

    end subroutine init_cc_amplitudes

    subroutine init_cc_hash(hash_table, hash_table_size) 
        ! routine to setup the hash for the quadrupels 
        type(cc_hash), pointer, intent(inout) :: hash_table(:)
        integer, intent(in) :: hash_table_size
        character(*), parameter :: this_routine = "init_cc_hash"

        integer :: i
        ! estimate the hash table size with num_doubles^2 for now.. 
        ! although this number can be reduced substantially i guess 

        allocate(hash_table(hash_table_size))

        do i = 1, hash_table_size
            hash_table(i)%found = .false. 
!             allocate(hash_table%ind(0,nifdbo))
!             hash_table(i)%ind = 0_n_int 
            hash_table(i)%amp = 0.0_dp 
            nullify(hash_table(i)%next)
        end do

    end subroutine init_cc_hash

    subroutine calc_n_triples() 
        ! this routine calculates the number of "important" triples 
        ! from the samples singles and doubles amplitudes: 
        ! essentiall calulating T1 * T2 
        character(*), parameter :: this_routine = "calc_n_triples"

        integer :: i, j, ia(2,1), jk_cd(2,2)

        do i = 1, cc_ops(1)%n_ops 
            ! for each t_i^a i have to check if a double excitation is 
            ! possible on top of it.. 
            if (.not. cc_ops(1)%set_flag(i)) cycle

            ia = cc_ops(1)%get_ex(i)

            do j = 1, cc_ops(2)%n_ops

                ! i also have to check if the amplitute was set 
                if (.not. cc_ops(2)%set_flag(j)) cycle 

                jk_cd = cc_ops(2)%get_ex(j)

                ! check if the operators fit..
                if (any(ia(1,1) == jk_cd(1,:)) .or. any(ia(2,1) == jk_cd(2,:))) then
                    cycle
                end if
                ! if it fits increase the triples counter:
                est_triples = est_triples + 1

            end do
        end do

    end subroutine calc_n_triples

    subroutine calc_n_quads 
        ! this routines calculates the number of "important" quadrupels
        ! calculating T2*T2 essentially 
        ! todo: do i double count here? and how do i solve this issue? 
        ! since different excitation generators can lead to the same 
        ! quadruple excitation or? 
        ! yes there are worst case 18 different permutations leading to the 
        ! same quadruple excitation! 
        ! so either divide by this number! 
        ! or, for not too big numbers of double excitations create a 
        ! specific hash table for the quadruple indexing out of the 
        ! 8 indices i,j,k,l and a,b,c,d .. 
        ! which would also allow me to store the estimated quadrupels 
        ! magnitute(when ignoring the influence of the T1 and T3 operators!)
        ! but when i store the tripels also i could also estimate the 
        ! quadrupels "exactly" .. which would be quite nice actually! 
        ! todo: hash-table or dividing by 18.. 
        character(*), parameter :: this_routine = "calc_n_quads" 

        integer :: i, j, ij_ab(2,2), kl_cd(2,2), ijab_klcd(8), hash_ind
        integer(n_int) :: temp_int(0:nifdbo)
        real(dp) :: phase, amp
        logical :: t_found

        if (t_store_hash_quadrupels) then 
          do i = 1, cc_ops(2)%n_ops 
                if (cc_ops(2)%set_flag(i)) then 

                    ij_ab = cc_ops(2)%get_ex(i)

                    amp = cc_ops(2)%get_amp(i)

                    do j = i + 1, cc_ops(2)%n_ops 

                        if (cc_ops(2)%set_flag(j)) then 

                            kl_cd = cc_ops(2)%get_ex(j)

                            ! can i just: any(a==b)
                            if (unique_quad_ind(ij_ab, kl_cd)) then

                                ! i also want to store the amplitudes 
                                ! i gues.. 
                                amp = amp * cc_ops(2)%get_amp(j)

                                ! here i have to first order the indices 
                                ! so i < j < k < l and a < b < c < d 
                                ! and find the phase factor between them 
                                call order_quad_indices(ij_ab, kl_cd, phase, &
                                    ijab_klcd)

                                amp = amp * phase 

                                ! and i need a unique quantity associated with 
                                ! these 8 orbital -> just set all the 
                                ! corresponding orbitals 
                                temp_int = ilutref(0:nifdbo,1)
                                clr_orb(temp_int, ij_ab(1,1))
                                clr_orb(temp_int, ij_ab(1,2))
                                clr_orb(temp_int, kl_cd(1,1))
                                clr_orb(temp_int, kl_cd(1,2))

                                set_orb(temp_int, ij_ab(2,1))
                                set_orb(temp_int, ij_ab(2,2))
                                set_orb(temp_int, kl_cd(2,1))
                                set_orb(temp_int, kl_cd(2,2))

                                call cc_hash_look_up(ijab_klcd, temp_int, quad_hash, &
                                    hash_ind, t_found)

                                ! it it is found in the hash-table i want to 
                                ! update the amplitude 
                                if (t_found) then 
                                    call cc_hash_update(quad_hash, hash_ind, & 
                                        temp_int, amp)
                                else 
                                    ! otherwise i want to add the entry 
                                    call cc_hash_add(quad_hash, hash_ind, & 
                                        temp_int, amp)

                                    ! and only in this case we found a new 
                                    ! quadruple coefficient
                                    est_quads = est_quads + 1
                                end if
                            end if
                        end if
                    end do
                end if
            end do

        else 
            do i = 1, cc_ops(2)%n_ops 
                if (cc_ops(2)%set_flag(i)) then 

                    ij_ab = cc_ops(2)%get_ex(i)

                    do j = i + 1, cc_ops(2)%n_ops 

                        if (cc_ops(2)%set_flag(j)) then 

                            kl_cd = cc_ops(2)%get_ex(j)

                            ! can i just: any(a==b)
                            if (unique_quad_ind(ij_ab, kl_cd)) then

                                est_quads = est_quads + 1

                            end if
                        end if
                    end do
                end if
            end do
        end if

    end subroutine calc_n_quads 

    subroutine cc_hash_look_up(ind, tgt, hash_table, hash_val, t_found)
        integer, intent(in) :: ind(:) 
        integer(n_int), intent(in) :: tgt(0:nifdbo)
        type(cc_hash), pointer :: hash_table(:)
        integer, intent(out) :: hash_val 
        logical, intent(out) :: t_found
        character(*), parameter :: this_routine = "cc_hash_look_up"

        type(cc_hash), pointer :: temp_node 

        t_found = .false. 

        hash_val = FindWalkerHash(ind, size(hash_table))

        temp_node => hash_table(hash_val) 

        if (temp_node%found) then 
            do while (associated(temp_node)) 
                if (all(tgt == temp_node%ind)) then 
                    t_found = .true. 
                    exit 
                end if
                temp_node => temp_node%next
            end do
        end if

        nullify(temp_node)

    end subroutine cc_hash_look_up

    subroutine cc_hash_update(hash_table, hash_val, tgt, amp)
        ! routine to update the amplitude of a found entry 
        type(cc_hash), pointer, intent(inout) :: hash_table(:)
        integer, intent(in) :: hash_val
        integer(n_int), intent(in) :: tgt(:) 
        real(dp), intent(in) :: amp 
        character(*), parameter :: this_routine = "cc_hash_update"

        type(cc_hash), pointer :: temp_node 
        logical :: found 

        found = .false.

        temp_node => hash_table(hash_val) 

        do while (associated(temp_node)) 
            if (all(temp_node%ind == tgt)) then 
                ! this is the correct entry! 
                found = .true.
                temp_node%amp = temp_node%amp + amp 
                exit 
            end if
            temp_node => temp_node%next
        end do

        ASSERT(found)

    end subroutine cc_hash_update

    subroutine cc_hash_add(hash_table, hash_val, tgt, amp) 
        ! this is a routine to add a hash table entry 
        ! it means the entry was not there yet or has to be added to 
        ! the linked list 
        type(cc_hash), pointer, intent(inout) :: hash_table(:)
        integer, intent(in) :: hash_val 
        integer(n_int), intent(in) :: tgt(:)
        real(dp), intent(in) :: amp 
        character(*), parameter :: this_routine = "cc_hash_add" 

        type(cc_hash), pointer :: temp_node 

        temp_node => hash_table(hash_val) 

        if (.not.temp_node%found) then 
            ! this means this entry is empty so we cann fill it here 
            temp_node%found = .true.
            allocate(temp_node%ind(0:nifdbo))
            temp_node%ind = tgt 
            temp_node%amp = amp 
        else 
            ! loop through the linked list 
            do while (associated(temp_node%next))
                temp_node => temp_node%next 
            end do
            allocate(temp_node%next) 
            nullify(temp_node%next%next)
            temp_node%next%found = .true.
            allocate(temp_node%next%ind(0:nifdbo))
            temp_node%next%ind = tgt 
            temp_node%next%amp = amp 

            ! for testing cound the number of clashes: 
            n_clashes = n_clashes + 1
        end if

        nullify(temp_node)

    end subroutine cc_hash_add


    subroutine order_quad_indices(ij_ab, kl_cd, phase, ijab_klcd)
        integer, intent(inout) :: ij_ab(2,2), kl_cd(2,2)
        real(dp), intent(out) :: phase 
        integer, intent(out) :: ijab_klcd(8)
        character(*), parameter :: this_routine = "order_quad_indices"

        integer :: ij(2), ab(2), kl(2), cd(2), i, j, k, l, a, b, c, d
        integer :: n 
        ! the correctly order t_ij^ab * t_kl^cd has the sign convention -1!
        ! and we can be sure that the inputed (ij),(ab) etc. are ordered i < j
        ! within the T2 operators, since only those get stored by convention 
        ! and also that i < k and a < c 

        ij = ij_ab(1,:)
        ab = ij_ab(2,:)
        kl = kl_cd(1,:)
        cd = kl_cd(2,:) 

        n = 1 

        ! assert i < k and a < c
        ASSERT(ij_ab(1,1) < kl_cd(1,1))
        ASSERT(ij(1) < ij(2))
        ASSERT(cd(1) < cd(2))
        ASSERT(ab(1) < ab(2))
        ASSERT(cd(1) < cd(2))

        ! with this: i am not so sure: a < c .. 
        ! ok i can not be sure that a < c! so i have to sort more here! 

        ! sort the electrons: 
        ! if j < k everything is fine and in order! 
        if (ij_ab(1,2) < kl_cd(1,1)) then 
            ! everything is fine.. 

        else if (ij_ab(1,2) < kl_cd(1,2)) then 
            ! if j < l then i have to switch j and k 
            call swap(ij(2), kl(1))
            n = n + 1 

        else 
            ! this means j > l so we have to swap j -> l and then k -> l 
            call swap(ij(2), kl(2))
            call swap(ij(2), kl(1))
            ! wich leaves the phase unchanged 

        end if 

        ! and check if we did everything correctly: 
        ASSERT(ij(1) < ij(2))
        ASSERT(kl(1) < kl(2))
        ASSERT(ij(2) < kl(1))

        ! i am not so sure any more if a < c anymore 
        ! now sort the orbitals but here a > c is also possible! 
        if (ij_ab(2,1) < kl_cd(2,1)) then 
            ! so a < c 
            if (ij_ab(2,2) < kl_cd(2,1)) then 
                ! b < c so in this case nothing has to be done! 


            else if (ij_ab(2,2) < kl_cd(2,2)) then 
                ! b < d: so b and c have to be swapped! 
                call swap(ab(2),cd(1))
                n = n + 1 

            else 
                ! b > d: so switch b -> d and d -> c 
                call swap(ab(2),cd(2))
                call swap(ab(2),cd(1))
                ! and this does not change the phase 
            end if
        else 
            ! so this means a > c which implies b > c also! 
            ! so if b < d we know everything 
            if (ij_ab(2,2) < kl_cd(2,2)) then 
                ! then c < a < b < d 
                call swap(ab(1),cd(1))
                call swap(ab(2),cd(1))
                ! so no phase factor! 

            else if (ij_ab(2,1) > kl_cd(2,2)) then 
                ! here b and a > d 
                ! so c < d < a < b 
                call swap(ab(1),cd(1))
                call swap(ab(2),cd(2))
                ! also here no phase factor! 

            else 
                ! this means c < a < d < b 
                call swap(ab(1),ab(2))
                ! ba, cd
                call swap(cd(1),cd(2))
                ! ba, dc 
                call swap(ab(1),cd(2))
                ! cadb
                ! here we get a phase 
                n = n + 1
            end if
        end if

        ASSERT(ab(1) < ab(2))
        ASSERT(cd(1) < cd(2))
        ASSERT(ab(2) < cd(1))

        ! now switch the indices
        ij_ab(1,:) = ij 
        ij_ab(2,:) = ab 
        kl_cd(1,:) = kl
        kl_cd(2,:) = cd 

        ! and also store a linear index 
        ijab_klcd = [ij,ab,kl,cd]

        ! and calculate the appropriate phase
        phase = (-1.0_dp)**n

    end subroutine order_quad_indices

    logical function unique_quad_ind(ij_ab, kl_cd)
        integer, intent(in) :: ij_ab(2,2), kl_cd(2,2)

        unique_quad_ind = all(ij_ab(1,1) /= kl_cd(1,:) .and. ij_ab(1,2) /= kl_cd(1,:) & 
                    .and. ij_ab(2,1) /= kl_cd(2,:) .and. ij_ab(2,2) /= kl_cd(2,:))

    end function unique_quad_ind

    ! do it other way.. only store the possible non-zero cluster operators! 
    subroutine fill_cc_amplitudes() 
        ! design decisions: since the singles are not so many in general 
        ! and because i could need them to correct the doubles amplitudes
        ! store all of the possible ones! and encode them specifically through 
        ! (i) and (a) 
        character(*), parameter :: this_routine = "fill_cc_amplitudes"
        
        integer :: idet, ic, ex(2,cc_order), j, i
        integer :: ia, ib, ja, jb, ind
        integer, allocatable :: n_excits(:)

        integer :: a, b, elec_i, elec_j, orb_a, orb_b 
        logical :: t_par, t_store_full_doubles = .true.

        real(dp) :: amp

        real(dp), allocatable :: temp_amps(:)
        integer, allocatable :: temp_ops(:,:,:) 

        integer :: ab(2), ac(2), bc(2), ij(2), ik(2), jk(2) 
        integer :: c, k, ij_ac, ij_bc, ij_ab, ik_bc, ik_ac, ik_ab
        integer :: jk_bc, jk_ac, jk_ab, elec_k, orb_c, jc, ka, kb, kc, n
        integer(n_int) :: temp_ilut(0:niftot)
        integer :: dummy_ind, dummy_hash, temp_nI(nel)
        logical :: tSuccess

        HElement_t(dp) :: sign_tmp(lenof_sign) 

        ! for this it is helpful to have an upper limit of the number of 
        ! possible amplitudes, but just do it for the singles for now..
        allocate(n_excits(2)) 

        n_excits = calc_number_of_excitations(nOccAlpha, nOccBeta, cc_order, & 
            nbasis/2)

        allocate(cc_ops(cc_order))

        ! and do a nice initialization depending on the order 
        do i = 1, cc_order 

            cc_ops(i)%order = i 

        end do

        allocate(cc_ops(1)%amplitudes(n_excits(1)))
        allocate(cc_ops(1)%operators(n_excits(1),2,1))
        allocate(cc_ops(1)%set_flag(n_excits(1)))

        cc_ops(1)%operators = 0 
        cc_ops(1)%amplitudes = 0.0_dp
        cc_ops(1)%set_flag = .false.

        cc_ops(1)%n_ops = n_excits(1)

        ! first figure out the number of double and fill the singles 
        ! amplitudes! 

        n_doubles = 0
        do idet = 1, int(totwalkers) 
            ! for now also count the exact number of triples and quadrupels
            ic = FindBitExcitLevel(ilutRef(:,1), CurrentDets(:,idet))
            select case (ic) 
            case (1) 
                n_singles = n_singles + 1 
                call extract_sign(CurrentDets(:,idet), sign_tmp)
                call get_bit_excitmat(iLutRef(:,1), CurrentDets(:,idet), ex, ic)
                ind = cc_ops(1)%get_ind(ex(1,1),ex(2,1))

                ! for now only do it for single runs.. do i need the normalising?
                amp = sign_tmp(1) / AllNoatHf(1)
                if (abs(amp) > EPS) then 
                    cc_ops(1)%amplitudes(ind) = amp
                    cc_ops(1)%operators(ind,:,1) = ex(:,1)
                    cc_ops(1)%set_flag(ind) = .true. 
                end if

            case (2) 
                n_doubles = n_doubles + 1

            case (3)
                n_triples = n_triples + 1

            case (4) 
                n_quads = n_quads + 1

            end select
        end do

        print *, "direct sampled doubles: ", n_doubles

        if (allocated(cc_ops(2)%operators))     deallocate(cc_ops(2)%operators)
        if (allocated(cc_ops(2)%amplitudes))    deallocate(cc_ops(2)%amplitudes)

        ! here i have a choice to only store the non-zero contributions 
        ! or encode them in the full-list to better access them later on. 
        ! especially when we want to calculate the coupled-cluster 
        ! triples and quadrupels.. 

        if (t_store_full_doubles) then 
            allocate(cc_ops(2)%operators(n_excits(2), 2, 2))
            allocate(cc_ops(2)%amplitudes(n_excits(2)))
            allocate(cc_ops(2)%set_flag(n_excits(2)))
            cc_ops(2)%n_ops = n_excits(2)
        else
            allocate(cc_ops(2)%operators(n_doubles, 2, 2))
            allocate(cc_ops(2)%amplitudes(n_doubles))
            allocate(cc_ops(2)%set_flag(n_doubles))
            cc_ops(2)%n_ops = n_doubles
        end if

        cc_ops(2)%operators = 0
        cc_ops(2)%amplitudes = 0.0_dp
        cc_ops(2)%set_flag = .false.
        
        do idet = 1, int(totwalkers)

            ic = FindBitExcitLevel(ilutRef(:,1), CurrentDets(:,idet))

            if (ic == 2) then 
                call extract_sign(CurrentDets(:,idet), sign_tmp)
                ! i hope this routine is not buggy..
                call get_bit_excitmat(iLutRef(:,1), CurrentDets(:,idet), ex, ic)

                ia = cc_ops(1)%get_ind(ex(1,1),ex(2,1))
                ib = cc_ops(1)%get_ind(ex(1,1),ex(2,2))
                ja = cc_ops(1)%get_ind(ex(1,2),ex(2,1))
                jb = cc_ops(1)%get_ind(ex(1,2),ex(2,2))

                ! check first if the amp gets 0 due to the singles
                amp = sign_tmp(1)/AllNoatHf(1) + & 
                    cc_ops(1)%get_amp(ja)*cc_ops(1)%get_amp(ib) -  & 
                    cc_ops(1)%get_amp(ia)*cc_ops(1)%get_amp(jb)

                if (abs(amp) > EPS) then 
                    if (t_store_full_doubles) then 
                        ind = cc_ops(2)%get_ind(ex(1,:),ex(2,:))
                    else 
                        ind = j 
                        j = j + 1
                    end if

                    cc_ops(2)%operators(ind,:,:) = ex(:,1:2)
                    cc_ops(2)%amplitudes(ind) = amp
                    cc_ops(2)%set_flag(ind) = .true.

                end if
            end if
        end do

        ! i just realise that i have to run over all the possible double 
        ! excitations, since t_ij^ab = c_ij^ab - t_i^a... could be 
        ! non-zero, just from the single contribution.. 
        if (t_store_full_doubles) then 
            do i = 1, nel 
                elec_i = projedet(i,1)
                do j = i + 1, nel 
                    elec_j = projedet(j,1)
                    t_par = same_spin(elec_i, elec_j) 

                    if (t_par) then 
                        do a = 1, nvirt 
                            orb_a = mask_virt_ni(a,1)
                            do b = a + 1, nvirt 
                                orb_b = mask_virt_ni(b,1) 

                                if (same_spin(orb_a, orb_b) .and. &
                                    same_spin(orb_a, elec_i)) then 

                                    ind = cc_ops(2)%get_ind([elec_i,elec_j], &
                                        [orb_a,orb_b]) 

                                    if (.not. cc_ops(2)%set_flag(ind)) then 

                                        ia = cc_ops(1)%get_ind([elec_i],[orb_a])
                                        ib = cc_ops(1)%get_ind([elec_i],[orb_b])
                                        ja = cc_ops(1)%get_ind([elec_j],[orb_a])
                                        jb = cc_ops(1)%get_ind([elec_j],[orb_b])

                                        
                                        amp = cc_ops(1)%get_amp(ja)*cc_ops(1)%get_amp(ib) -  & 
                                              cc_ops(1)%get_amp(ia)*cc_ops(1)%get_amp(jb)

                                        if (abs(amp) > EPS) then 
                                            cc_ops(2)%operators(ind,1,:) = [elec_i, elec_j]
                                            cc_ops(2)%operators(ind,2,:) = [orb_a, orb_b] 
                                            cc_ops(2)%amplitudes(ind) = amp 
                                            cc_ops(2)%set_flag(ind) = .true. 

                                            n_doubles = n_doubles + 1
                                        end if
                                    end if
                                end if
                            end do
                        end do
                    else
                        do a = 1, nvirt 
                            orb_a = mask_virt_ni(a,1)
                            do b = a + 1, nvirt 
                                orb_b = mask_virt_ni(b,1) 

                                if (.not. same_spin(orb_a, orb_b)) then

                                    ind = cc_ops(2)%get_ind([elec_i,elec_j], &
                                        [orb_a,orb_b]) 

                                    if (.not. cc_ops(2)%set_flag(ind)) then 

                                        ia = cc_ops(1)%get_ind([elec_i],[orb_a])
                                        ib = cc_ops(1)%get_ind([elec_i],[orb_b])
                                        ja = cc_ops(1)%get_ind([elec_j],[orb_a])
                                        jb = cc_ops(1)%get_ind([elec_j],[orb_b])

                                        
                                        amp = cc_ops(1)%get_amp(ja)*cc_ops(1)%get_amp(ib) -  & 
                                              cc_ops(1)%get_amp(ia)*cc_ops(1)%get_amp(jb)

                                        if (abs(amp) > EPS) then 
                                            cc_ops(2)%operators(ind,1,:) = [elec_i, elec_j]
                                            cc_ops(2)%operators(ind,2,:) = [orb_a, orb_b] 
                                            cc_ops(2)%amplitudes(ind) = amp 
                                            cc_ops(2)%set_flag(ind) = .true. 

                                            n_doubles = n_doubles + 1
                                        end if
                                    end if
                                end if
                            end do
                        end do
                    end if
                end do
            end do
        end if
        
        print *, "doubles after singles contribution: ", n_doubles

        ! and should i also do the triples.. just to check maybe? 
        ! but here i definetly only want to store the non-zero 
        ! contributions.. but.. i actually would need to run over 
        ! all the possible ones also since the contribution could be 
        ! approximated by T1*T2 or T1^3.. 

        ! i could loop.. but i do not want to encode all of them! 
        ! but essentially i have to first allocate a vector of all possible 
        ! ones to be sure.. 
        if (cc_order == 3) then 
            allocate(temp_amps(n_excits(3))) 
            allocate(temp_ops(n_excits(3),2,3))

            n = 1

            ! first analyze the wavefunction: 
            do idet = 1, int(totwalkers)
               ic = FindBitExcitLevel(iLutRef(:,1), CurrentDets(:,idet))

               if (ic == 3) then 

                   call extract_sign(CurrentDets(:,idet), sign_tmp)
                   call get_bit_excitmat(iLutRef(:,1), CurrentDets(:,idet), ex, ic)

                   ia = cc_ops(1)%get_ind(ex(1,1),ex(2,1))
                   ib = cc_ops(1)%get_ind(ex(1,1),ex(2,2))
                   ic = cc_ops(1)%get_ind(ex(1,1),ex(2,3))

                   ja = cc_ops(1)%get_ind(ex(1,2),ex(2,1))
                   jb = cc_ops(1)%get_ind(ex(1,2),ex(2,2))
                   jc = cc_ops(1)%get_ind(ex(1,2),ex(2,3))

                   ka = cc_ops(1)%get_ind(ex(1,3),ex(2,1))
                   kb = cc_ops(1)%get_ind(ex(1,3),ex(2,2))
                   kc = cc_ops(1)%get_ind(ex(1,3),ex(2,3))

                   ij = [ex(1,1),ex(1,2)]
                   ik = [ex(1,1),ex(1,3)]
                   jk = [ex(1,2),ex(1,3)]

                   ab = [ex(2,1),ex(2,2)]
                   ac = [ex(2,1),ex(2,3)]
                   bc = [ex(2,2),ex(2,3)]

                   jk_bc = cc_ops(2)%get_ind(jk,bc)
                   jk_ac = cc_ops(2)%get_ind(jk,ac)
                   jk_ab = cc_ops(2)%get_ind(jk,ab)

                   ik_bc = cc_ops(2)%get_ind(ik,bc)
                   ik_ac = cc_ops(2)%get_ind(ik,ac)
                   ik_ab = cc_ops(2)%get_ind(ik,ab)

                   ij_bc = cc_ops(2)%get_ind(ij,bc)
                   ij_ac = cc_ops(2)%get_ind(ij,ac)
                   ij_ab = cc_ops(2)%get_ind(ij,ab)

                   amp = sign_tmp(1)/AllNoatHf(1) &
                       - cc_ops(1)%get_amp(ia)*cc_ops(2)%get_amp(jk_bc) &
                       + cc_ops(1)%get_amp(ib)*cc_ops(2)%get_amp(jk_ac) & 
                       - cc_ops(1)%get_amp(ic)*cc_ops(2)%get_amp(jk_ab) & 
                       + cc_ops(1)%get_amp(ja)*cc_ops(2)%get_amp(ik_bc) & 
                       - cc_ops(1)%get_amp(jb)*cc_ops(2)%get_amp(ik_ac) & 
                       + cc_ops(1)%get_amp(jc)*cc_ops(2)%get_amp(ik_ab) & 
                       - cc_ops(1)%get_amp(ka)*cc_ops(2)%get_amp(ij_bc) & 
                       + cc_ops(1)%get_amp(kb)*cc_ops(2)%get_amp(ij_ac) & 
                       - cc_ops(1)%get_amp(kc)*cc_ops(2)%get_amp(ij_ab) & 
                       - cc_ops(1)%get_amp(ia)*cc_ops(1)%get_amp(jb)*cc_ops(1)%get_amp(kc) & 
                       + cc_ops(1)%get_amp(ia)*cc_ops(1)%get_amp(jc)*cc_ops(1)%get_amp(kb) & 
                       + cc_ops(1)%get_amp(ib)*cc_ops(1)%get_amp(ja)*cc_ops(1)%get_amp(kc) & 
                       - cc_ops(1)%get_amp(ib)*cc_ops(1)%get_amp(jc)*cc_ops(1)%get_amp(ka) & 
                       - cc_ops(1)%get_amp(ic)*cc_ops(1)%get_amp(ja)*cc_ops(1)%get_amp(kb) & 
                       + cc_ops(1)%get_amp(ic)*cc_ops(1)%get_amp(jb)*cc_ops(1)%get_amp(ka)

                   if (abs(amp) > EPS) then 

                       temp_amps(n) = amp 
                       temp_ops(n,:,:) = ex

                       n = n + 1

                   end if
               end if
            end do
            print *, "direct sampled triples: ", n - 1

            do i = 1, nel 
                elec_i = projedet(i,1)
                do j = i + 1, nel 
                    elec_j = projedet(j,1)
                    ij = [elec_i, elec_j]
                    do k = j + 1, nel 
                        elec_k = projedet(k,1)

                        ik = [elec_i,elec_k] 
                        jk = [elec_j,elec_k]

                        do a = 1, nvirt
                            orb_a = mask_virt_ni(a,1)

                            ia = cc_ops(1)%get_ind([elec_i],[orb_a])
                            ja = cc_ops(1)%get_ind([elec_j],[orb_a])
                            ka = cc_ops(1)%get_ind([elec_k],[orb_a])

                            do b = a + 1, nvirt 
                                orb_b = mask_virt_ni(b,1)
                                ab = [orb_a,orb_b] 
                                ib = cc_ops(1)%get_ind([elec_i],[orb_b])
                                jb = cc_ops(1)%get_ind([elec_j],[orb_b])
                                kb = cc_ops(1)%get_ind([elec_k],[orb_b])

                                ij_ab = cc_ops(2)%get_ind(ij,ab)
                                ik_ab = cc_ops(2)%get_ind(ik,ab) 
                                jk_ab = cc_ops(2)%get_ind(jk,ab)

                                do c = b + 1, nvirt 
                                    orb_c = mask_virt_ni(c,1)
                                    ac = [orb_a,orb_c]
                                    bc = [orb_b,orb_c]

                                    ic = cc_ops(1)%get_ind([elec_i],[orb_c])
                                    jc = cc_ops(1)%get_ind([elec_j],[orb_c])
                                    kc = cc_ops(1)%get_ind([elec_k],[orb_c]) 

                                    ij_ac = cc_ops(2)%get_ind(ij,ac)
                                    ik_ac = cc_ops(2)%get_ind(ik,ac)
                                    jk_ac = cc_ops(2)%get_ind(jk,ac)

                                    ij_bc = cc_ops(2)%get_ind(ij,bc)
                                    ik_bc = cc_ops(2)%get_ind(ik,bc)
                                    jk_bc = cc_ops(2)%get_ind(jk,bc)

                                    ! and now i have to check if the 
                                    ! triple excitation (ijk|abc) is in the 
                                    ! wallker list! use the hash-table!
                                    ! todo 
                                    temp_ilut = iLutRef(:,1) 
                                    ! make the excitation and search! 
                                    clr_orb(temp_ilut, elec_i)
                                    clr_orb(temp_ilut, elec_j)
                                    clr_orb(temp_ilut, elec_k)
                                    set_orb(temp_ilut, orb_a)
                                    set_orb(temp_ilut, orb_b)
                                    set_orb(temp_ilut, orb_c)

                                    call decode_bit_det(temp_nI, temp_ilut)

                                    call hash_table_lookup(temp_nI, temp_ilut, &
                                        nifdbo, HashIndex, CurrentDets, dummy_ind, &
                                        dummy_hash, tSuccess)

                                    if (.not. tSuccess) then
                                        amp = -cc_ops(1)%get_amp(ia)*cc_ops(2)%get_amp(jk_bc) &
                                           + cc_ops(1)%get_amp(ib)*cc_ops(2)%get_amp(jk_ac) & 
                                           - cc_ops(1)%get_amp(ic)*cc_ops(2)%get_amp(jk_ab) & 
                                           + cc_ops(1)%get_amp(ja)*cc_ops(2)%get_amp(ik_bc) & 
                                           - cc_ops(1)%get_amp(jb)*cc_ops(2)%get_amp(ik_ac) & 
                                           + cc_ops(1)%get_amp(jc)*cc_ops(2)%get_amp(ik_ab) & 
                                           - cc_ops(1)%get_amp(ka)*cc_ops(2)%get_amp(ij_bc) & 
                                           + cc_ops(1)%get_amp(kb)*cc_ops(2)%get_amp(ij_ac) & 
                                           - cc_ops(1)%get_amp(kc)*cc_ops(2)%get_amp(ij_ab) & 
                                           - cc_ops(1)%get_amp(ia)*cc_ops(1)%get_amp(jb)*cc_ops(1)%get_amp(kc) & 
                                           + cc_ops(1)%get_amp(ia)*cc_ops(1)%get_amp(jc)*cc_ops(1)%get_amp(kb) & 
                                           + cc_ops(1)%get_amp(ib)*cc_ops(1)%get_amp(ja)*cc_ops(1)%get_amp(kc) & 
                                           - cc_ops(1)%get_amp(ib)*cc_ops(1)%get_amp(jc)*cc_ops(1)%get_amp(ka) & 
                                           - cc_ops(1)%get_amp(ic)*cc_ops(1)%get_amp(ja)*cc_ops(1)%get_amp(kb) & 
                                           + cc_ops(1)%get_amp(ic)*cc_ops(1)%get_amp(jb)*cc_ops(1)%get_amp(ka)

                                       if (abs(amp) > EPS) then 
                                           ! then add this to the list! 
                                           temp_amps(n) = amp
                                           temp_ops(n,1,:) = [elec_i,elec_j,elec_k]
                                           temp_ops(n,2,:) = [orb_a,orb_b,orb_c]
                                           n = n + 1
                                       end if
                                   end if
                                end do
                            end do
                        end do
                    end do
                end do
            end do

            print *, "triples after T1*T2 and T1^3: ", n - 1
            ! and then allocate the actual array: 
            allocate(cc_ops(3)%amplitudes(n-1))
            allocate(cc_ops(3)%operators(n-1,2,3))

            cc_ops(3)%amplitudes = temp_amps(1:n-1)
            cc_ops(3)%operators = temp_ops(1:n-1,:,:)
            cc_ops(3)%n_ops = n - 1
        end if


    end subroutine fill_cc_amplitudes

    function get_ex(this, ind) result(ex) 
        ! this function gives the specific excitation, if present in the 
        ! cc_amplitudes, which are encoded in a linear fashion 
        ! with the convention that all the electron and orbital indices are 
        ! always provided in an ordered(from lowest to highest) fashion
        class(cc_amplitude), intent(in) :: this
        integer, intent(in) :: ind
        integer :: ex(2, this%order)
        character(*), parameter :: this_routine = "get_ex"

        integer :: ij, ab, cum, i, j, a, b, nij

        ASSERT(ind > 0) 
        select case (this%order) 
        case (1) 

            ex = this%operators(ind,:,:)
            return

            ! fix this below if i want to actually encode that directly and 
            ! not store the operators..
            ASSERT(ind <= nel * (nbasis - nel))

            ! it is a single excition so this is easy to decompose 
            ! this decomposition depends on the encoding below!
            ex(2,1) = mod(ind-1,(nbasis - nel))
            ! lets hope the integer division is done correctly.. on all compilers
            ex(1,1) = int((ind - ex(2,1))/(nbasis - nel)) + 1
            ! and modify by nel the orbital index again to get the real 
            ! orbital index! 
            ex(2,1) = ex(2,1) + nel + 1

        case (2) 

            ex = this%operators(ind,:,:)
            return

            ! fix this below if i want to actually encode that directly and 
            ! not store the operators..
            call stop_all(this_routine, "still buggy for double excitations!")
            ! first we have to get ij, ab back: 

            nij = nel * (nel - 1) / 2

            ab = mod(ind - 1, nij)
            ij = (ind - ab)/(nij) + 1 

            ab = ab + 1 

            ! then we have to get ij and ab from those indices in the same 
            ! way 
            ! but here it gets more tricky because we can just do the mod.. 
            ! maybe i have to store the indices in the end.. 
            i = 1
            cum = (nel - i)

            print *, ""
            print *, "nij: ", nij 
            print *, "ab: ", ab 
            print *, "ij: ", ij 
            print *, "cum: ", cum

            ! do a stupid sum..
            do while(cum < ij .and. i <= nel)
                i = i + 1
                cum = cum + (nel - i)
            end do

            j = ij + i - (cum - nel + i)

            ! and the same for ab 
            a = 1 
            cum = (nbasis - nel - a)
            do while (cum < ab .and. a <= (nbasis - nel)) 
                a = a + 1
                cum = cum + (nbasis - nel - a)
            end do

            b = ab + a - (cum - (nbasis - nel) + a)

            ex(1,:) = [i,j]
            ex(2,:) = [a,b] + nel

            print *, "(i,j): ", i,j 
            print *, "(a,b): ", a, b
        case default 
            call stop_all(this_routine, "higher than cc_order 2 not yet implemented!")

        end select

    end function get_ex

    function get_ind(this, elec_ind, orb_ind) result(ind)
        ! depending on the cc_order this encodes all the electron and 
        ! orbital indices of an excitation in a unique linear fashion 
        ! we can assume the electron and orbital indices to be ordered from 
        ! highest to lowest
        class(cc_amplitude), intent(in) :: this 
        integer, intent(in) :: elec_ind(this%order), orb_ind(this%order) 
        integer :: ind

        character(*), parameter :: this_routine = "get_ind"

        integer :: ij, ab, nij, orb(2)

        ind = -1 

        ! to i want to access this with invalid excitations?
        ASSERT(all(elec_ind > 0))
!         ASSERT(all(elec_ind <= nel))
!         ASSERT(all(orb_ind > nel))
        ASSERT(all(orb_ind <= nbasis))

        ! and assert ordered inputs.. 
        ! or do we want to order it here, if it is not ordered? 
        ! tbd!
        ASSERT(minloc(elec_ind,1) == 1) 
        ASSERT(maxloc(elec_ind,1) == this%order)

        ASSERT(minloc(orb_ind,1) == 1) 
        ASSERT(maxloc(orb_ind,1) == this%order)

        select case (this%order)
        case (1) 
            ! single excitation 
            ! the elec_ind goes from 1 to nel and the orb_ind goes from 
            ! nel + 1 to nbasis (has nbasis - nel) possible values 
            ! can we assume a closed shell ordered reference? 
            ! with the nel lowest orbitals occupied? 
            ! otherwise this is a bit more tricky.. 
            ind = ind_matrix_singles(elec_ind(1), orb_ind(1))
            return

            ind = (elec_ind(1) - 1) * (nbasis - nel) + (orb_ind(1) - nel)

        case (2) 

            ij = linear_elec_ind(elec_ind(1), elec_ind(2))
            ab = linear_orb_ind(orb_ind(1), orb_ind(2))

            if (ij == 0 .or. ab == 0) then 
                ind = 0
            else
                ind = ind_matrix_doubles(ij,ab)
            end if

            return

            call stop_all(this_routine, "still buggy for double excitations!")
            ! double excitation 
            ! first encode the ij electron index in a linear fashion
            ASSERT(elec_ind(1) < elec_ind(2))
            ASSERT(orb_ind(1) < orb_ind(2))

            orb = orb_ind - nel 

            ij = (elec_ind(1) - 1) * nel - elec_ind(1)*(elec_ind(1) - 1)/2  + (elec_ind(2) - elec_ind(1))
!             ij = (elec_ind(2) - 1)*elec_ind(2)/2 + elec_ind(1) - 1
            ! i shift the orb_indices by nel.. 
!             ab = (orb_ind(1) - nel - 1)*(nbasis - nel) + (orb_ind(2) - orb_ind(1))
            ab = (orb(1) - 1) *(nbasis - nel) - orb(1)*(orb(1) - 1)/2 + (orb(2) - orb(1))
!             ab = (orb_ind(2) - nel - 1)*(orb_ind(2)-nel)/2 + orb_ind(1) - nel -1

            nij = nel * (nel - 1) / 2 
            ind = (ij - 1) * nij + ab 

        case default 
            call stop_all(this_routine, "higher than cc_order 2 not yet implemented!")
        end select

    end function get_ind

!     subroutine calc_cc_amplitudes
!         character(*), parameter :: this_routine = "calc_cc_amplitudes" 
! 
!         integer :: idet, i
!         integer, allocatable :: n_excits(:)
!         type(cc_amplitude), allocatable :: cc_amp(:)
! 
!         ! i want to calculate the amplitudes up to a certain order given 
!         ! in the input 
! 
!         ! for this it is helpful to have an upper limit of the number of 
!         ! possible amplitudes 
!         allocate(n_excits(cc_order)) 
! 
!         n_excits = calc_number_of_excitations(nOccAlpha, nOccBeta, cc_order, & 
!             nbasis/2)
! 
!         allocate(cc_amp(cc_order))
! 
!     end subroutine calc_cc_amplitudes


    subroutine dongxia_amplitudes 
 
! dongxia test
        integer :: IC,idet,ierr,Nvirt,ind,idi,idj,idk,ida,idb,idc
        integer :: Ex(2,3)
        real(dp), allocatable :: C1(:),C2(:)
        real(dp), allocatable :: T1(:),T2(:)
        real(dp) :: sign_tmp(lenof_sign),C3,T3

          Ex = 0

! Dongxia 1.8.2017
! to obtain CC amplitudes from walker numbers
        Nvirt = Nbasis - Nel
        allocate(C1(Nvirt*nel),stat=ierr)
          c1 = 0.0_dp
        allocate(C2(Nvirt*(Nvirt+1)*nel*(nel+1)/4),stat=ierr)
          c2 = 0.0_dp
        allocate(T1(Nvirt*nel),stat=ierr)
          t1 = 0.0_dp
        allocate(T2(Nvirt*(Nvirt+1)*nel*(nel+1)/4),stat=ierr)
          t2 = 0.0_dp
        write(6,*) 'dongxia testing C1 and T1'
! look for C1 and obtain T1
        do idet=1, int(totwalkers)
           IC = 4
           call get_bit_excitmat(iLutRef,CurrentDets(:,idet),ex,IC)
           if (IC==1) then
              call extract_sign(CurrentDets(:,idet),sign_tmp)
              idi = Ex(1,1)
              ida = Ex(2,1) - Nel
              ind = ind1(idi,ida)
              C1(ind) = sign_tmp(1) / AllNoatHf(1)
              T1(ind) = C1(ind)
              if(abs(c1(ind)) > 1.0d-3) &
               write(6,'(3I5,F20.12)') ex(1,1),ex(2,1),ind, C1(ind)
           end if
        end do
! look for C2, and obtain T2
        write(6,*)'dongxia testing double excitations'
        do idet=1, int(totwalkers)
           IC = 4
           call get_bit_excitmat(iLutRef,CurrentDets(:,idet),ex,IC)
           if (IC==2) then
              call extract_sign(CurrentDets(:,idet),sign_tmp)
              idi = Ex(1,1)
              ida = Ex(2,1)-Nel
              idj = Ex(1,2)
              idb = Ex(2,2)-Nel
              ind = ind2(idi,idj,ida,idb)
              C2(ind) = sign_tmp(1)/AllNoatHF(1)
              T2(ind) = C2(ind)+T1(ind1(idj,ida))*T1(ind1(idi,idb)) &
                        -T1(ind1(idi,ida))*T1(ind1(idj,idb))
              if(abs(c2(ind)) > 1.0d-3) &
               write(6,'(4I10,2F20.12)') idi,idj,ida,idb,C2(ind),T2(ind)
           end if
        end do
! look for C3, and compare it with T3 (from T1^3, T1T2)
        write(6,*)'dongxia testing triple excitations'
        do idet = 1, int(totwalkers)
           IC = 4
           call get_bit_excitmat(iLutRef,CurrentDets(:,idet),ex,IC)
           if(IC==3)then
             c3 = 0.0_dp
             t3 = 0.0_dp
             call extract_sign(CurrentDets(:,idet),sign_tmp)
             idi = Ex(1,1)
             idj = Ex(1,2)
             idk = Ex(1,3)
             ida = Ex(2,1) - Nel
             idb = Ex(2,2) - Nel
             idc = Ex(2,3) - Nel
             C3 = sign_tmp(1)/AllNoatHF(1)
! calculate t3 (by t1^3, t1^t2 and t2^t1)
! contribution from t1^3/6.
             t3 = t3 &
                 +t1(ind1(idi,ida))*t1(ind1(idj,idb))*t1(ind1(idk,idc)) &
                 +t1(ind1(idi,idb))*t1(ind1(idj,idc))*t1(ind1(idk,ida)) &
                 +t1(ind1(idi,idc))*t1(ind1(idj,ida))*t1(ind1(idk,idb)) &
                 -t1(ind1(idi,ida))*t1(ind1(idj,idc))*t1(ind1(idk,idb)) &
                 -t1(ind1(idi,idc))*t1(ind1(idj,idb))*t1(ind1(idk,ida)) &
                 -t1(ind1(idi,idb))*t1(ind1(idj,ida))*t1(ind1(idk,idc)) 
! contribution from t1^t2 (and t2^t1, which is the same)
             t3 = t3 &
                 +t1(ind1(idi,ida))*t2(ind2(idj,idk,idb,idc)) &
                 -t1(ind1(idi,idb))*t2(ind2(idj,idk,ida,idc)) &
                 +t1(ind1(idi,idc))*t2(ind2(idj,idk,ida,idb)) &
                 -t1(ind1(idj,ida))*t2(ind2(idi,idk,idb,idc)) &
                 +t1(ind1(idj,idb))*t2(ind2(idi,idk,ida,idc)) &
                 -t1(ind1(idj,idc))*t2(ind2(idi,idk,ida,idb)) &
                 +t1(ind1(idk,ida))*t2(ind2(idi,idj,idb,idc)) &
                 -t1(ind1(idk,idb))*t2(ind2(idi,idj,ida,idc)) &
                 +t1(ind1(idk,idc))*t2(ind2(idi,idj,ida,idb))  

                      
             if(abs(c3)+abs(t3) > 1.0d-3) &
              write(6,'(3I5,4X,3I5,3F20.12)')idi,idj,idk,ida,idb,idc,c3,t3,t3/c3-1
           end if
        end do

        deallocate(T2,stat=ierr)
        deallocate(T1,stat=ierr)
        deallocate(C2,stat=ierr)
        deallocate(C1,stat=ierr)
! end
    end subroutine dongxia_amplitudes

      pure integer function ind1(i,a)

      implicit none
      integer, intent(in) :: i, a

      ind1 = (i-1)*(nbasis-nel)+a
      
      return
      end function ind1

      pure integer function ind2(i,j,a,b)

      implicit none
      integer, intent(in) :: i,j,a,b
      integer :: ij, ab, nvirt,nab
      
      ij = (j-1)*j/2 + i
      ab = (b-1)*b/2 + a
      nvirt = nbasis - nel
      nab = nvirt*(nvirt+1)/2
       
      ind2 = (ij-1)*nab + ab

      return
      end function ind2      

    function calc_number_of_excitations(n_alpha, n_beta, max_excit, n_orbs) &
            result(n_excits)
        ! i might want a routine which calculates the number of all possible 
        ! excitations for each excitation level 
        integer, intent(in) :: n_alpha, n_beta, max_excit, n_orbs
        integer :: n_excits(max_excit)

        integer :: n_parallel(max_excit,2), i, j, k

        ! the number of all possible single excitations, per spin is: 
        
        n_parallel = 0
        n_excits = 0

        do i = 1, max_excit
            n_parallel(i,1) = calc_n_parallel_excitations(n_alpha, n_orbs, i)
            n_parallel(i,2) = calc_n_parallel_excitations(n_beta, n_orbs, i)
        end do

        ! i can set this up in a general way.. 

        do i = 1, max_excit

            ! the 2 parallel spin species always are added 
            n_excits(i) = sum(n_parallel(i,:)) 

            ! and i have to zig-zag loop over the other possible combinations 
            ! to lead to this level of excitation.. 
            j = i - 1 
            k = 1 

            do while (j >= k) 

                n_excits(i) = n_excits(i) + n_parallel(j,1) * n_parallel(k,2) 

                ! and if j /= k do it the other way around too 
                if (j /= k) then 
                    n_excits(i) = n_excits(i) + n_parallel(j,2) * n_parallel(k,1)

                end if
                ! and in/decrease the counters 
                j = j - 1 
                k = k + 1

                ! in this way i calculate eg.: 
                ! n_ij^ab = n_ij_ab(u + d) + n_i^a(u) * n_i^a(d) 
                ! and 
                ! n_ijk^abc = n_ijk^abc(u + d) + n_ij^ab(u) * n_i^a(d) and spin-flipped..

            end do
        end do

    end function calc_number_of_excitations

    function calc_n_single_excits(n_elecs, n_orbs) result(n_single_excits)
        ! this calculates the number of possible single excitations for a 
        ! given spin species of electrons and the number of spatial orbitals!
        integer, intent(in) :: n_elecs, n_orbs
        integer :: n_single_excits

        ! it is the number of electrons which can be excited and the number 
        ! of available orbitals for this spin!
!         n_single_excits = n_elecs * (n_orbs - n_elecs)

        ! with binomial it is just: 
        n_single_excits = binomial(n_elecs, 1) * binomial(n_orbs - n_elecs, 1)

    end function calc_n_single_excits

    function calc_n_parallel_excitations(n_elecs, n_orbs, ic) result(n_parallel)
        ! this function determines the number of parallel spin excitaiton 
        ! with each electrons having the same spin for a given excitation 
        ! level and number of electrons and available orbitals for this spin 
        integer, intent(in) :: n_elecs, n_orbs, ic 
        integer :: n_parallel

        n_parallel = binomial(n_elecs, ic) * binomial(n_orbs - n_elecs, ic) 

    end function calc_n_parallel_excitations

    integer function binomial(n, k)
        ! write a new binomial function using the fortran2008 standard 
        ! gamma function 
        integer, intent(in) :: n, k 
#ifdef __DEBUG
        character(*), parameter :: this_routine = "binomial"
#endif

        ASSERT(n >= 0)
        ASSERT(k >= 0)

        if (k > n) then 
            binomial = 0
        else
            binomial = nint(gamma(real(n) + 1) / (gamma(real(n - k) + 1) * gamma(real(k) + 1)))
        end if

    end function binomial 

end module cc_amplitudes
