#include "macros.h"

! i guess for a new module it would be nice to have one module which tests 
! all the functions and programs of a module .. check that out.. 

! all the global variables in NECI make it a bit tough to do unit-tests

! no.. i guess, atleast in my framework it would be better to write one 
! unit test module which tests all the functionality of a other module 
program test_k_space_hubbard 
    
    use k_space_hubbard

    use constants 

    use fruit 

    use SystemData, only: t_k_space_hubbard, t_lattice_model, nel, nbasis, & 
                          t_trans_corr, G1, nBasisMax, nOccBeta, nOccAlpha, & 
                          bhub, uhub, omega, trans_corr_param_2body, & 
                          t_trans_corr, t_trans_corr_2body, trans_corr_param, & 
                          thub, tpbc, treal, ttilt, TSPINPOLAR, & 
                          tCPMD, tVASP, tExch, tHphf, tNoSymGenRandExcits, tKPntSym, &
                          t_twisted_bc, twisted_bc, arr, brr

    use bit_rep_data, only: niftot, nifd

    use lattice_mod, only: lat, lattice, get_helement_lattice_general, & 
                           get_helement_lattice_ex_mat, get_helement_lattice

    use dsfmt_interface, only: dsfmt_init

    use OneEInts, only: GetTMatEl, tOneElecDiag, tCPMDSymTMat

    use procedure_pointers, only: get_umat_el

    use IntegralsData, only: umat

    use DetBitOps, only: EncodeBitDet

    use fcimcdata, only: pDoubles, pParallel

    use sort_mod, only: sort

    use util_mod, only: choose, get_free_unit

    use bit_reps, only: decode_bit_det, encode_sign

    use SymExcitDataMod, only: kTotal

    use lanczos_wrapper, only: frsblk_wrapper

    use unit_test_helpers

    use lattice_models_utils, only: gen_all_excits_k_space_hubbard, &
                                    gen_all_triples_k_space, create_hilbert_space_kspace, &
                                    gen_all_doubles_k_space

    use analyse_wf_symmetry, only: analyze_full_wavefunction_sym, &
                                   t_symmetry_mirror, t_symmetry_rotation

    use bit_reps, only: init_bit_rep

    implicit none 

    integer :: failed_count 


    call init_fruit()
    call dsfmt_init(0)

    ! misuse the unit tests for now to also do an exact study.. 
    call exact_study() 
    ! run the test-driver 
    call k_space_hubbard_test_driver()
    call fruit_summary()
    call fruit_finalize() 

    call get_failed_count(failed_count) 
    if (failed_count /= 0) stop -1 

contains 

    subroutine exact_study

        use DetCalcData, only: nkry, nblk, b2l, ncycle
        integer, allocatable :: hilbert_space(:,:), nI(:), nJ(:)
        integer(n_int), allocatable :: hilbert_space_ilut(:,:)
        real(dp) :: tmp_sign(lenof_sign)
        real(dp) :: J, U
        real(dp), allocatable :: J_vec(:), U_vec(:)
        integer :: i, n_eig, m
        real(dp), allocatable :: e_values(:), e_vecs(:,:)
        real(dp) :: k1(2),k2(2), k_vec(3)
        integer ::  n_excits, k, iunit, l, iunit2
        integer(n_int), allocatable :: excits(:,:)
        real(dp) :: sum_doubles, sum_doubles_trans
        character(30) :: filename, U_str, nel_str, lattice_name, J_str
        logical :: t_do_diags, t_do_subspace_study, t_input_U, t_U_vec
        HElement_t(dp), allocatable :: hamil(:,:)
        integer, allocatable :: sort_ind(:)
        real(dp), allocatable :: twist_x_vec(:), twist_y_vec(:)
        logical :: t_do_twisted_bc, t_twisted_vec, t_ignore_k, t_do_ed
        real(dp), allocatable :: epsilon_kvec(:), left_ev(:,:)
        integer, allocatable :: ind(:)
        real(dp) :: E_ref, gap
        logical :: t_analyse_gap, t_J_vec, t_input_J, t_sym
        real(dp), allocatable :: sym_labels(:), sym_labels_2(:)

        tmp_sign = 0.0_dp

        t_do_diags = .false.
        t_do_subspace_study = .true.
        t_input_U = .false.
        t_input_J = .false.
        t_U_vec = .true.
        t_J_vec = .true.
        t_do_twisted_bc = .false.
        t_twisted_vec = .false.
        t_analyse_gap = .false.
        t_ignore_k = .false.
        t_do_ed = .true.

        call init_k_space_unit_tests()

        ! i have to define the lattice here.. 
        lat => lattice('tilted', 3, 3, 1,.true.,.true.,.true.,'k-space')

!         x = [(-lat%dispersion_rel_orb(i), i = 1, 24)]
!         ind = [(i, i = 1, 24)]
! 
!         call sort(x,ind) 
! 
!         k1 = [2*pi/8., 10.0*pi/(3.0*8.0)]
!         k2 = [-2*pi/8.0,2*pi/8.0]
!         print *, "k | ole k |  e(k): "
!         do i = 1, 24 
!             k_vec = lat%get_k_vec(ind(i))
!             print *,lat%get_k_vec(ind(i)),"|", k_vec(1)*k1 + k_vec(2)*k2, "|", x(i)
!         end do

        if (t_input_U) then
            print *, "input U: "
            read(*,*) U
        else if (t_U_vec) then
            U_vec = linspace(1.0,8.0,2)
        else
            U = 4.0
        end if

        if (t_input_J) then
            print *, "input J:"
            read(*,*), J
        else if (t_J_vec) then
            J_vec = linspace(0.0,1.0,2)
        else 
            J = 0.1
        end if
        
        nel = 14
        allocate(nI(nel))
        allocate(nJ(nel))
        nj = 0

        nI = [(i, i = 1,nel)]
        ! 46 in 50:
!         nI = [ 9,10,11,12,13,14,15,16,21,22,23,24,25,26,27,28,29,30,&
!             37,38,39,40,41,42,43,44,45,46,55,56,57,58,59,60,61,62,63,64,73,74,75,76,77,78,79,80]

        ! 14 in 18:
!         nI = [ 3,4, 5, 6, 11,12,13,14,15,16,23,24,25,26 ]

        ! 6 in 9 k = 1 1
!         nI = [3,4,7,8,9,10]
        ! 6 in 9 k = 0 0
!         nI = [3,4,9,10,15,16]
        ! 10 in 9, k = 0:
!         nI = [3,4,7,8,9,10,11,12,15,16]
        ! 10 in 9, k = 1,-1:
!         nI = [1,2,3,4,7,8,9,10,11,12]

        ! chain:
        ! 6 in 6, k = 0
!         nI = [3,4,5,6,7,8]

        ! 6 in 6, k = 3
!         nI = [2,3,4,5,6,7]

        ! 4 in 6, k = 0
!         nI = [3,5,6,8]

        ! 4 in 6, k = -2
!         nI = [3,4,5,6]

        ! 4 in 6, k = 2
!         nI = [5,6,7,8]

        ! tilted 2x2:
        ! 8 in 8, k = 0
!         ni = [1,2,3,4,5,6,7,8,11,12] 


!         ni = [7,8,9,10,11,12,13,14,19,20,21,22,23,24,25,26,33,34,35,36,37,38,39,40]
!         ni = [7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]
!         nI = [3,4,5,6,7,8,11,12,13,14,15,16,21,22,23,24,25,26]
!         nI = [(i, i = 1,nel)]
!         nI = [9,10,11,12,13,14,15,16,17,18,21,22,23,24,25,26,27,28,&
!             29,30,37,38,39,40,41,42,43,44,45,46,55,56,57,58,59,60,61,62,63,64,&
!             71,72,73,74,75,76,77,78,79,80]
!         ni = [7,8,15,16,17,18]
!         nI = [21,23,24,26]
!         nI = [15,16]

!         nI = [1,2,3,4,5,6]
!         nI = [5,6,7,8,9,10]

        ! setup lanczos:
        nblk = 4
        nkry = 8 
        ncycle = 200
        b2l = 1.0e-13_dp

        print *, "nkry: ", nkry
        print *, "nblk: ", nblk
        print *, "b2l: ", b2l
        print *, "ncycle: ", ncycle

        call init_bit_rep()

        ! use twisted bc in this case.. 
        if (t_do_twisted_bc) then
            t_twisted_bc = .true. 
            if (t_twisted_vec) then 
                twist_x_vec = linspace(0.0,1.0,1000)
!                 allocate(twist_x_vec(1))
!                 twist_x_vec = 0.0
!                 allocate(twist_y_vec(1))
!                 twist_y_vec = 0.0
                twist_y_vec = linspace(0.0,1.0,1000)


                if (t_analyse_gap) then

                    call setup_system(lat, nI, J, U)
                    iunit = get_free_unit()
                    open(iunit, file = 'gap_vs_k')
                    write(iunit,*) '# x-twist, y-twist, gap'

                    do i = 1, size(twist_x_vec)
                        twisted_bc(1) = twist_x_vec(i)
                        twisted_bc(2) = 0.0
!                         twisted_bc(2) = 2*twist_y_vec(i)
                        epsilon_kvec = [(-lat%dispersion_rel_orb(l), l = 1, nBasis/2)] 

                        call sort(epsilon_kvec)
                        gap = epsilon_kvec(nel/2 + 1) - epsilon_kvec(nel/2)

                        write(iunit,*) twisted_bc(1:2), gap

                    end do
                    close(iunit)
                    call stop_all("here", "now")
                end if

                call setup_system(lat, nI, J, U, hilbert_space)
                n_eig = 5
                iunit = get_free_unit()
                open(iunit, file = 'one_body_vs_twist')
                iunit2 = get_free_unit()
                open(iunit2, file = 'energy_vs_twist')
                write(iunit,*) '# x-twist, y-twist, E(k)'
                write(iunit2,*) '# x-twist, y-twist, k-total, E_ref, energies'
                allocate(e_values(n_eig))
                allocate(e_vecs(size(hilbert_space,2),n_eig))
                allocate(epsilon_kvec(nBasis/2))
                allocate(ind(nBasis/2))

                do i = 1, size(twist_x_vec)
                    do k = 1, size(twist_y_vec)
                        twisted_bc(1) = twist_x_vec(i)
                        twisted_bc(2) = twist_y_vec(k)

                        epsilon_kvec = [(-lat%dispersion_rel_orb(l), l = 1, nBasis/2)] 

                        write(iunit,*) twisted_bc(1:2), epsilon_kvec

                        if (t_ignore_k) then
                            ! if we want to fill up the nel lowest 
                            ! orbitals always.. 
                            ind = [(i, i = 1, nBasis/2)]
                            call sort(epsilon_kvec, ind)

                            nI = [2*ind(1:nel/2)-1,2*ind(1:nel/2)]

                            call sort(nI)

                            call setup_system(lat, nI, J, U, hilbert_space)

                            if (allocated(e_vecs)) deallocate(e_vecs)
                            allocate(e_vecs(size(hilbert_space,2),n_eig))

                        else 
                            call setup_system(lat, nI, J, U)
                        end if
                        E_ref = get_helement_lattice(nI,nI)
                        ! and also calculate the GS energy
                        ! do i need to setup the system more thoroughly for the 
                        ! GS energy? i think so..

                        print *, "k-twist: ", twisted_bc(1:2)
                        print *, "ni: ", nI
                        print *, "k-total: ", kTotal

                        if (t_do_ed) then
                            call frsblk_wrapper(hilbert_space, size(hilbert_space, 2), n_eig, e_values, e_vecs)
                        end if
                        write(iunit2,*) twisted_bc(1:2),kTotal, E_ref, e_values

                    end do
                end do
                close(iunit)
                close(iunit2)
                call stop_all("here", "now")
            else
                twisted_bc = 0.5
            end if
        end if
            


!         NIfTot = 0
!         nifd = 0

        if (t_do_diags) then
            t_trans_corr_2body = .true.
            trans_corr_param_2body = J_vec((1))
            three_body_prefac = 2.0_dp * (cosh(trans_corr_param_2body) - 1.0_dp) / real(omega**2,dp)
            call setup_system(lat, nI, J, U)

    !         call gen_all_doubles_k_space(nI, n_excits, excits)
            call gen_all_excits_k_space_hubbard(nI, n_excits, excits)

            write(U_str, *) int(U)
            write(nel_str, *) nel 
            lattice_name = "18_"

            filename = 'H_elements_' // trim(adjustl(nel_str)) // 'in' // &
                trim(adjustl(lattice_name)) // 'U_' // trim(adjustl(U_str))

            iunit = get_free_unit()

            open(iunit, file = filename)

            print *, "writing to file..."
            write(iunit, *) "# J, H diag, <I|H|K>, <K|H|I>:"
            do i = 1, size(J_vec)
                trans_corr_param_2body = J_vec((i))
                three_body_prefac = 2.0_dp * (cosh(trans_corr_param_2body) - 1.0_dp) / real(omega**2,dp)
     
                sum_doubles = 0.0_dp
                sum_doubles_trans = 0.0_dp

                do k = 1, n_excits
                    call decode_bit_det(nJ, excits(:,k))

                    sum_doubles = sum_doubles + get_helement_lattice(nI, nJ)
                    sum_doubles_trans = sum_doubles_trans + get_helement_lattice(nJ,nI)

                end do

                write(iunit, *)  J_vec(i), get_diag_helement_k_sp_hub(nI), sum_doubles, sum_doubles_trans
            end do

            close(iunit)
            t_trans_corr_2body = .false.

            call stop_all("here", "now")
        end if

        if (t_do_subspace_study) then
            allocate(hilbert_space(nel,8))
            n_eig = 8

            if (lat%get_nsites() == 50) then
                ! 50-site:
                ! closed shell:
                ! 9 10 79 80
                hilbert_space(:,1) = [ 9,10,11,12,13,14,15,16,21,22,23,24,25,26,27,28,29,30,&
                    37,38,39,40,41,42,43,44,45,46,55,56,57,58,59,60,61,62,63,64,73,74,75,76,77,78,79,80]

                ! 17 18 71 72
                hilbert_space(:,2) = [ 11,12,13,14,15,16,17,18,21,22,23,24,25,26,27,28,29,30,&
                    37,38,39,40,41,42,43,44,45,46,55,56,57,58,59,60,61,62,63,64,71,72,73,74,75,76,77,78]


                ! open shell:
                ! 9 17 72 80
                hilbert_space(:,3) = [ 9,11,12,13,14,15,16,17,21,22,23,24,25,26,27,28,29,30,&
                    37,38,39,40,41,42,43,44,45,46,55,56,57,58,59,60,61,62,63,64,72,73,74,75,76,77,78,80]

                ! 9 18 71 80
                hilbert_space(:,4) = [9 ,11,12,13,14,15,16,18,21,22,23,24,25,26,27,28,29,30,&
                    37,38,39,40,41,42,43,44,45,46,55,56,57,58,59,60,61,62,63,64,71,73,74,75,76,77,78,80]

                ! 9 18 72 79
                hilbert_space(:,5) = [ 9,11,12,13,14,15,16,18,21,22,23,24,25,26,27,28,29,30,&
                    37,38,39,40,41,42,43,44,45,46,55,56,57,58,59,60,61,62,63,64,72,73,74,75,76,77,78,79]

                ! 10 17 71 80
                hilbert_space(:,6) = [ 10,11,12,13,14,15,16,17,21,22,23,24,25,26,27,28,29,30,&
                    37,38,39,40,41,42,43,44,45,46,55,56,57,58,59,60,61,62,63,64,71,73,74,75,76,77,78,80]

                ! 10 17 72 79
                hilbert_space(:,7) = [ 10,11,12,13,14,15,16,17,21,22,23,24,25,26,27,28,29,30,&
                    37,38,39,40,41,42,43,44,45,46,55,56,57,58,59,60,61,62,63,64,72,73,74,75,76,77,78,79]
                
                ! 10 18 71 79
                hilbert_space(:,8) = [10 ,11,12,13,14,15,16,18,21,22,23,24,25,26,27,28,29,30,&
                    37,38,39,40,41,42,43,44,45,46,55,56,57,58,59,60,61,62,63,64,71,73,74,75,76,77,78,79]

            else if (lat%get_nsites() == 18) then

                ! 18-site 
                ! closed shell:
                ! 3 4 25 26
                hilbert_space(:,1) = [ 3,4, 5, 6, 11,12,13,14,15,16, 23,24,25,26 ]
                ! 7 8 21 22 
                hilbert_space(:,2) = [  5, 6,7,8 , 11,12,13,14,15,16,21,22 ,23,24 ]
                ! open shell: 
                ! 3 7 22 26
                hilbert_space(:,3) = [ 3, 5, 6, 7, 11,12,13,14,15,16,22 ,23,24,26 ]
                ! 3 8 21 26
                hilbert_space(:,4) = [ 3, 5, 6, 8, 11,12,13,14,15,16,21 ,23,24,26 ]
                ! 3 8 22 25
                hilbert_space(:,5) = [ 3, 5, 6,8 , 11,12,13,14,15,16,22 ,23,24,25 ]
                ! 4 7 21 26
                hilbert_space(:,6) = [ 4, 5, 6,7 , 11,12,13,14,15,16,21 ,23,24,26 ]
                ! 4 7 22 25
                hilbert_space(:,7) = [ 4, 5, 6,7 , 11,12,13,14,15,16, 22,23,24,25 ]
                ! 4 8 21 25
                hilbert_space(:,8) = [ 4, 5, 6, 8, 11,12,13,14,15,16,21 ,23,24,25 ]
            end if

            allocate(hilbert_space_ilut(0:niftot,size(hilbert_space,2)))

            hilbert_space_ilut = 0_n_int

            do i = 1, size(hilbert_space,2)
                call EncodeBitDet(hilbert_space(:,i), hilbert_space_ilut(:,i))
            end do

            t_symmetry_rotation = .true.
            t_symmetry_mirror = .true.
            allocate(e_values(n_eig))
            allocate(e_vecs(n_eig,n_eig))
            allocate(left_ev(n_eig,n_eig))

            allocate(sort_ind(n_eig))
            t_trans_corr_2body = .true. 

            iunit = get_free_unit()
            open(iunit, file = 'basis')
            call print_matrix(transpose(hilbert_space), iunit)
            close(iunit)

            if (t_U_vec)then
                do i = 1, size(U_vec)
                    U = U_vec(i)
                    write(U_str,*) int(U)

                    if (t_J_vec) then
                        do k = 1, size(J_vec)
                            J = J_vec(k)

                            write(J_str,*) J
                            filename = "energy_e_vectors_U_" // trim(adjustl(U_str)) & 
                                // "_J_" // trim(adjustl(J_str))

                            iunit = get_free_unit()

                            call setup_system(lat, nI, J, U)

                            hamil = create_hamiltonian(hilbert_space)

                            t_sym = check_symmetric(hamil)
                            if (t_sym) then
                                call eig_sym(hamil, e_values, e_vecs) 
                            else
                                call eig(hamil, e_values, e_vecs) 
                            end if

                            print *, "basis:" 
                            call print_matrix(transpose(hilbert_space))
                            print *, "hamil: "
                            call print_matrix(hamil)

                            sort_ind = [(i, i = 1,n_eig)]
                            call sort(e_values, sort_ind)
                            print *, "e-values: ", e_values
                            print *, "sort_ind: ", sort_ind
                            e_vecs = e_vecs(:,sort_ind)
                            print *, "e-vecs: "
                            call print_matrix(e_vecs)
                            if (.not. t_sym) then
                                call eig(hamil, e_values, left_ev, .true.) 
                                sort_ind = [(i, i = 1,n_eig)]
                                call sort(e_values, sort_ind)

                                left_ev = left_ev(:,sort_ind)
                            end if

                            ! also do the symmetry study.. 

                            open(iunit, file = filename)

                            if (t_sym) then
                                write(iunit,*) &
                                    "# U, J, E, EV, 1, 90, 180, 270, m_x, m_y, m_d, m_o, S"
                            else
                                write(iunit,*) &
                                    "# U, J, E, left_ev, right_ev, 0, 90, 180, 270, m_x, m_y, m_d, m_o, S"
                            end if

                            do l = 1, n_eig
                                do m = 1, size(hilbert_space,2)
                                    tmp_sign(1) = e_vecs(m,l)
                                    call encode_sign(hilbert_space_ilut(:,m), tmp_sign)
                                end do
                                call analyze_full_wavefunction_sym(sym_labels, hilbert_space_ilut)
                                if (.not. t_sym) then
                                    do m = 1, size(hilbert_space,2)
                                        tmp_sign(1) = left_ev(m,l)
                                        call encode_sign(hilbert_space_ilut(:,m), tmp_sign)
                                    end do
                                    call analyze_full_wavefunction_sym(sym_labels_2, hilbert_space_ilut)
                                end if

                                if (t_sym) then
                                    write(iunit,*) U, J, e_values(l), e_vecs(:,l), nint(sym_labels)
                                else
                                    write(iunit,*) &
                                        U, J, e_values(l), e_vecs(:,l), left_ev(:,l), nint(sym_labels)
                                end if
                            end do
                            close(iunit)

                        end do
                    end if
                end do

                call stop_all("here", "now")

            else

                call setup_system(lat, nI, J, U)

                hamil = create_hamiltonian(hilbert_space)

                call eig(hamil, e_values, e_vecs) 

                print *, "basis:" 
                call print_matrix(transpose(hilbert_space))
                print *, "hamil: "
                call print_matrix(hamil)

                sort_ind = [(i, i = 1,8)]
                call sort(e_values, sort_ind)
                print *, "e-values: ", e_values
                print *, "sort_ind: ", sort_ind
                e_vecs = e_vecs(:,sort_ind)
                print *, "e-vecs: "
                call print_matrix(e_vecs)
            end if

        else
            call setup_system(lat, nI, J, U, hilbert_space)
            n_eig = 5

!             hamil = create_hamiltonian(hilbert_space)
!             allocate(e_values(n_eig))
!             allocate(e_vecs(n_eig, size(hilbert_space,2)))
!             call eig(hamil, e_values, e_vecs) 
!             sort_ind = [(i, i = 1,n_eig)]
!             ind = minloc(e_values,1) 
!             call sort(e_values, sort_ind)
!             e_vecs = e_vecs(sort_ind,:)

        end if

        print *, "k-vector : ", kTotal
! 
        allocate(e_values(n_eig))
        allocate(e_vecs(size(hilbert_space,2),n_eig))

        print *, "size hilbert: ", size(hilbert_space, 2)

        if (t_do_ed) then
            if (t_U_vec) then 
                U = U_vec(1)
                ! setup hilbert once
                iunit = get_free_unit()
                open(iunit, file = 'energy_vs_U')
                call setup_system(lat, nI, J, U, hilbert_space)
                write(iunit,*) "# U, E_0, E_1, E_2, E_3, E_4"
                do i = 1, size(U_vec)
                    print *, "U: ", U_vec(i)
                    uhub  = U_vec(i)
                    umat = h_cast(real(uhub,dp)/real(omega,dp))
                    call frsblk_wrapper(hilbert_space, size(hilbert_space, 2), n_eig, e_values, e_vecs)
                    write(iunit,*) U_vec(i), e_values
                end do
                close(iunit)
                call stop_all("here", "now")
            end if

            ! try too big systems here: 
            call frsblk_wrapper(hilbert_space, size(hilbert_space, 2), n_eig, e_values, e_vecs)
! 
            print *, "e_value lanczos:", e_values
            print *, "|i>, c_i:"
            do i = 1, size(hilbert_space,2)
                print *, hilbert_space(:,i), e_vecs(i,1), e_vecs(i,2), e_vecs(i,3)
            end do
        end if
        call stop_all("here", "now")
        call exact_transcorrelation(lat, nI, [J], U, hilbert_space) 

        call stop_all("here", "now")

    end subroutine exact_study
    
    subroutine setup_system(in_lat, nI, J, U, hilbert_space) 
        class(lattice), intent(in) :: in_lat
        integer, intent(in) :: nI(:) 
        real(dp), intent(in) :: J, U
        integer, intent(out), allocatable, optional :: hilbert_space(:,:)

        integer :: i, n_states
        integer(n_int), allocatable :: dummy(:,:) 

        bhub = -1.0_dp
        nel = size(nI)

        nOccBeta = 0 
        nOccAlpha = 0 

        do i = 1, nel 
            if (is_beta(nI(i)))  nOccBeta = nOccBeta + 1
            if (is_alpha(nI(i))) nOccAlpha = nOccAlpha + 1
        end do

        call setup_all(in_lat, J, U) 

        call setup_k_total(nI) 

        call init_dispersion_rel_cache()
        get_helement_lattice_ex_mat => get_helement_k_space_hub_ex_mat
        get_helement_lattice_general => get_helement_k_space_hub_general
        call init_tmat_kspace(lat)


        if (present(hilbert_space)) then
!             call create_hilbert_space(nI, n_states, hilbert_space, dummy, gen_all_excits_k_space_hubbard) 
            ! change to my new hilbert space creator 
            call create_hilbert_space_kspace(nOccAlpha, nOccBeta, in_lat%get_nsites(), & 
                nI, n_states, hilbert_space, dummy)
        end if

    end subroutine setup_system

    subroutine init_k_space_unit_tests()
        ! since there is so much annoying outside dependency, mainly due to 
        ! momentum-conservation symmetry, lets just initialize all 
        ! necessary stuff here in front so testing is not so annoying.. 

        print *, ""
        print *, "initializing k_space unit tests"
        print *, "for simplicity do the unit tests on a 4 electron system" 
        print *, "in 4 spatial orbitals" 
        nel = 4 
        nbasis = 8 
        NIfTot = 0
        nifd = 0 

        ! set the appropriate flags: 
        thub = .true.
        tpbc = .true. 
        treal = .false. 
        ! todo: also test for tilted lattices! since Ali likes them so much.. 
        TSPINPOLAR = .false. 
        tOneElecDiag = .false.

        tKPntSym = .true.
        ttilt =  .false. 
        tCPMDSymTMat = .false.
        tvasp = .false.
        tCPMD = .false. 

        tExch = .true. 

        thphf = .false.

        bhub = -1.0

        tNoSymGenRandExcits = .true. 

        lat => lattice('chain', 4, 1, 1,.true.,.true.,.true.,'k-space')

        t_k_space_hubbard = .true. 
        t_lattice_model = .true. 

        ! setup nBasisMax and also the same for nBasisMax
        call setup_nbasismax(lat)

        call setup_arr_brr(lat) 

        ! setup G1 properly
        ! do a function like: which depending on the lattice sets up everything
        ! necessary for this type of lattice! yes!
        call setup_g1(lat) 

        ! also need the tmat ready.. 
        call init_tmat_kspace(lat)
!         call setup_tmat_k_space(lat)

!         call setup_kPointToBasisFn(lat)

        call setup_k_space_hub_sym(lat) 

        ! and i also have to setup the symmetry table... damn.. 
        ! i have to setup umat also or
        uhub = 1.0
        omega = 4.0

        ! and i have to allocate umat.. 
        allocate(umat(1))
        umat = h_cast(real(uhub,dp)/real(omega,dp))

        get_umat_el => get_umat_kspace

        trans_corr_param_2body = 0.1
        three_body_prefac = 2.0_dp * (cosh(trans_corr_param_2body) - 1.0_dp) / real(omega**2,dp)

        ! also initialize the lattice get_helement pointers to use 
        ! detham to do the Lanczos procedure for bigger systems.. 
        call init_dispersion_rel_cache()
        get_helement_lattice_ex_mat => get_helement_k_space_hub_ex_mat
        get_helement_lattice_general => get_helement_k_space_hub_general
        call init_tmat_kspace(lat)

    end subroutine init_k_space_unit_tests

    subroutine k_space_hubbard_test_driver() 
        ! with all the annying symmetry stuff to set up, testing the 
        ! k-space hubbard is really annoying.. 
        ! this is the main function which calls all the other tests 
       
        call run_test_case(setup_g1_test, "setup_g1_test")
        call run_test_case(setup_nbasismax_test, "setup_nbasismax_test")
        call run_test_case(setup_tmat_k_space_test, "setup_tmat_k_space_test")
        call run_test_case(setup_kPointToBasisFn_test, "setup_kPointToBasisFn_test")

        call init_k_space_unit_tests()
        call run_test_case(test_3e_4orbs_par, "test_3e_4orbs_par")
        call run_test_case(test_3e_4orbs_trip, "test_3e_4orbs_trip")
        call run_test_case(test_4e_ms1, "test_4e_ms1")
        call run_test_case(test_4e_ms0_mom_1, "test_4e_ms0_mom_1")
        call run_test_case(test_3e_ms1, "test_3e_ms1")
!         call run_test_case(test_8e_8orbs, "test_8e_8orbs")
        call run_test_case(test_general, "test_general")

        call run_test_case(get_diag_helement_k_sp_hub_test, "get_diag_helement_k_sp_hub_test")
        call run_test_case(get_offdiag_helement_k_sp_hub_test, "get_offdiag_helement_k_sp_hub_test")
        call run_test_case(get_helement_k_space_hub_test, "get_helement_k_space_hub_test")
        call run_test_case(pick_spin_opp_elecs_test, "pick_spin_opp_elecs_test")
        call run_test_case(pick_from_cum_list_test, "pick_from_cum_list_test")
        call run_test_case(create_ab_list_hubbard_test, "create_ab_list_hubbard_test")
        call run_test_case(pick_ab_orbitals_hubbard_test, "pick_ab_orbitals_hubbard_test")
        call run_test_case(calc_pgen_k_space_hubbard_test, "calc_pgen_k_space_hubbard_test")
        call run_test_case(gen_excit_k_space_hub_test, "gen_excit_k_space_hub_test")
        call run_test_case(pick_three_opp_elecs_test, "pick_three_opp_elecs_test")
        call run_test_case(pick_spin_par_elecs_test, "pick_spin_par_elecs_test")
        call run_test_case(pick_a_orbital_hubbard_test, "pick_a_orbital_hubbard_test")
        call run_test_case(pick_bc_orbitals_hubbard_test,"pick_bc_orbitals_hubbard_test")
        call run_test_case(create_ab_list_par_hubbard_test, "create_ab_list_par_hubbard_test")
        call run_test_case(pick_ab_orbitals_par_hubbard_test, "pick_ab_orbitals_par_hubbard_test")
        call run_test_case(get_transferred_momenta_test, "get_transferred_momenta_test")
        call run_test_case(create_bc_list_hubbard_test, "create_bc_list_hubbard_test")
        call run_test_case(get_3_body_helement_ks_hub_test, "get_3_body_helement_ks_hub_test")
        call run_test_case(check_momentum_sym_test, "check_momentum_sym_test")
        call run_test_case(find_minority_spin_test, "find_minority_spin_test")
        call run_test_case(calc_pgen_k_space_hubbard_transcorr_test, "calc_pgen_k_space_hubbard_transcorr_test")
        call run_test_case(calc_pgen_k_space_hubbard_par_test, "calc_pgen_k_space_hubbard_par_test")
        call run_test_case(calc_pgen_k_space_hubbard_triples_test, "calc_pgen_k_space_hubbard_triples_test")
        call run_test_case(make_triple_test, "make_triple_test")
        call run_test_case(make_double_test, "make_double_test")
        call run_test_case(three_body_transcorr_fac_test, "three_body_transcorr_fac_test")
        call run_test_case(two_body_transcorr_factor_test, "two_body_transcorr_factor_test")
        call run_test_case(epsilon_kvec_test, "epsilon_kvec_test")
        call run_test_case(same_spin_transcorr_factor_test, "same_spin_transcorr_factor_test")
        call run_test_case(get_one_body_diag_test, "get_one_body_diag_test")
        call run_test_case(gen_parallel_double_hubbard_test, "gen_parallel_double_hubbard_test")
        call run_test_case(gen_triple_hubbard_test, "gen_triple_hubbard_test")
        call run_test_case(gen_excit_k_space_hub_test_stochastic, "gen_excit_k_space_hub_test_stochastic") 
        call run_test_case(gen_excit_k_space_hub_transcorr_test_stoch, "gen_excit_k_space_hub_transcorr_test_stoch")

    end subroutine k_space_hubbard_test_driver

    subroutine test_3e_4orbs_par

        integer :: hilbert_nI(3,6), i, j, work(18), info, nI(3), n_states
        HElement_t(dp) :: hamil(6,6), hamil_trancorr(6,6), tmp_hamil(6,6)
        real(dp) :: ev_real(6), ev_cmpl(6), left_ev(1,6), right_ev(1,6)
        real(dp) :: t_mat(6,6), trans_mat(6,6)
        integer, allocatable :: test_hilbert(:,:)
        integer(n_int), allocatable :: dummy(:,:) 

        nOccBeta = 2 
        nOccAlpha = 1
        nel = 3

        hilbert_nI(:,1) = [1,4,5]!
        hilbert_nI(:,2) = [2,3,5]!
        hilbert_nI(:,3) = [1,3,6]!
        hilbert_nI(:,4) = [3,7,8]!
        hilbert_nI(:,5) = [1,2,7]!
        hilbert_nI(:,6) = [5,6,7]!

        nI = [1,4,5]
        call setup_k_total(nI)

        call create_hilbert_space(nI, n_states, test_hilbert, dummy, gen_all_excits_k_space_hubbard)

        print *, "n_states: ", n_states 
        call print_matrix(transpose(test_hilbert))
        ! nice.. it acutally works and gets all states.. 

        t_trans_corr_2body = .false. 
        print *, "un-transcorrelated Hamiltonian: "

        hamil = create_hamiltonian(hilbert_nI)

        call print_matrix(hamil)

        t_mat = get_tranformation_matrix(hamil,2)

        trans_mat = matmul(matmul(matrix_exponential(-t_mat),hamil),matrix_exponential(t_mat))
        
        ! use the lapack routines to solve these quickly.. 
        print *, "eigen-values: ", calc_eigenvalues(hamil)

        print *, "transcorrelated Hamiltonian: "
        t_trans_corr_2body = .true.
        
        hamil_trancorr = create_hamiltonian(hilbert_nI)

        call print_matrix(hamil_trancorr)

        print *, "eigen-values: ", calc_eigenvalues(hamil_trancorr)

        print *, "transformed hamiltonian" 
        call print_matrix(trans_mat)

        print *, "eigen-values: ", calc_eigenvalues(trans_mat)
        
    end subroutine test_3e_4orbs_par

    subroutine test_3e_4orbs_trip

        integer :: hilbert_nI(3,6), i, j, work(18), info
        HElement_t(dp) :: hamil(6,6), hamil_trancorr(6,6), tmp_hamil(6,6)
        real(dp) :: ev_real(6), ev_cmpl(6), left_ev(1,6), right_ev(1,6)
        real(dp) :: t_mat(6,6), trans_mat(6,6)

        nOccBeta = 1 
        nOccAlpha = 2
        nel = 3

        hilbert_nI(:,1) = [2,3,4]
        hilbert_nI(:,2) = [1,2,6]
        hilbert_nI(:,3) = [4,5,8]
        hilbert_nI(:,4) = [4,6,7]
        hilbert_nI(:,5) = [2,7,8]
        hilbert_nI(:,6) = [3,6,8]

        ! first create the non-transcorrelated Hamiltonian
        t_trans_corr_2body = .false. 
        print *, "un-transcorrelated Hamiltonian: "
        hamil = create_hamiltonian(hilbert_nI)

        call print_matrix(hamil)
        
        ! use the lapack routines to solve these quickly.. 
        print *, "eigen-values: ", calc_eigenvalues(hamil)

        print *, "transcorrelated Hamiltonian: "
        t_trans_corr_2body = .true.
        hamil_trancorr = create_hamiltonian(hilbert_nI)

        call print_matrix(hamil_trancorr)

        print *, "eigen-values: ", calc_eigenvalues(hamil_trancorr)

        t_mat = get_tranformation_matrix(hamil,2) 
        trans_mat = matmul(matmul(matrix_exponential(-t_mat),hamil),matrix_exponential(t_mat))

        print *, "transformed hamiltonian" 
        call print_matrix(trans_mat)

        print *, "eigen-values: ", calc_eigenvalues(trans_mat)
        
    end subroutine test_3e_4orbs_trip

    subroutine test_4e_ms0_mom_1

        integer :: hilbert_nI(4,8)
        HElement_t(dp) :: hamil(8,8),hamil_trancorr(8,8), tmp_hamil(8,8) 
        real(dp) :: ev_real(8), ev_cmpl(8), left_ev(1,8), right_ev(1,8) 
        integer :: work(24), info, n
        real(dp) :: t_mat(8,8), trans_mat(8,8) 

        nOccBeta = 2 
        nOccAlpha = 2
        nel = 4 

        hilbert_nI(:,1) = [3,4,5,8] 
        hilbert_nI(:,2) = [1,2,4,5] 
        hilbert_nI(:,3) = [3,4,6,7]
        hilbert_nI(:,4) = [1,4,7,8]
        hilbert_nI(:,5) = [1,5,6,8]
        hilbert_nI(:,6) = [2,3,7,8]
        hilbert_nI(:,7) = [1,2,3,6]
        hilbert_nI(:,8) = [2,5,6,7]

        t_trans_corr_2body = .false. 

        hamil = create_hamiltonian(hilbert_nI)

        t_mat = get_tranformation_matrix(hamil,4) 

        trans_mat = matmul(matmul(matrix_exponential(-t_mat),hamil),matrix_exponential(t_mat))

        t_trans_corr_2body = .true. 

        hamil_trancorr = create_hamiltonian(hilbert_nI) 

        print *, "un-correlated hamiltonian: "
        call print_matrix(hamil) 

        print *, "eigenvalues: ", calc_eigenvalues(hamil)

        print *, "neci-correlated hamiltonian: "
        call print_matrix(hamil_trancorr)
        print *, "eigenvalues: ", calc_eigenvalues(hamil_trancorr)

        print *, "transformed hamiltonian" 
        call print_matrix(trans_mat)
        print *, "eigenvalues: ", calc_eigenvalues(trans_mat)

    end subroutine test_4e_ms0_mom_1
    subroutine test_4e_ms1

        integer :: hilbert_nI(4,4), i, j, three_e(3,3) 
        integer(n_int) :: hilbert_ilut(0:niftot,4)
        HElement_t(dp) :: hamil(4,4), hamil_trancorr(4,4), tmp_hamil(4,4)
        HElement_t(dp) :: tmp_3(3,3), hamil_3(3,3), hamil_3_trans(3,3)
        real(dp) :: ev_real(4), ev_cmpl(4)
        real(dp) :: left_ev(1,4), right_ev(1,4)
        real(dp) :: shift_12, shift_34, two_body, three_body, three_body_1, three_body_2
        real(dp) :: two_body_1, two_body_2
        integer :: work(12), info
        real(dp) :: t_mat(4,4), trans_mat(4,4)

        nOccBeta = 3 
        nOccAlpha = 1 
        nel = 4

        hilbert_nI(:,1) = [1,2,3,5]
        hilbert_nI(:,2) = [3,4,5,7]
        hilbert_nI(:,3) = [1,3,7,8] 
        hilbert_nI(:,4) = [1,5,6,7]

        ! first create the non-transcorrelated Hamiltonian
        t_trans_corr_2body = .false. 
        print *, "un-transcorrelated Hamiltonian: "
        
        hamil = create_hamiltonian(hilbert_nI)
        call print_matrix(hamil)

        ! use the lapack routines to solve these quickly.. 
        print *, "eigen-values: ", calc_eigenvalues(hamil)

        print *, "transcorrelated Hamiltonian: "
        t_trans_corr_2body = .true.

        hamil_trancorr = create_hamiltonian(hilbert_nI)
        call print_matrix(hamil_trancorr)

        print *, "eigen-values: ", calc_eigenvalues(hamil_trancorr)

        t_mat = get_tranformation_matrix(hamil,3) 

        trans_mat = matmul(matmul(matrix_exponential(-t_mat),hamil),matrix_exponential(t_mat))

        print *, "transformed hamiltonian: "
        call print_matrix(trans_mat)

        print *, "eigen-values: ", calc_eigenvalues(trans_mat)

        ! i know now, where there is an error-- between those 2 matrix 
        ! elements: 
!         print *, "============== Excitation:  ======================"
!         print *, "nI:", hilbert_nI(:,1)
!         print *, "nJ:", hilbert_nI(:,2)
!         print *, "H_ij: ", hamil_trancorr(1,2), hamil_trancorr(2,1)
! 
!         ! the excitation is from (1,2) -> (4,7)
!         print *, "excitation: (1,2) -> (4,7)" 
!         ! ex: should be + 
!         ! 1 2 
!         ! 4 7 
!         ! total 1 2 3 5 should be - 
! 
!         two_body_1 = two_body_transcorr_factor(G1(1)%k, G1(7)%k)! + & 
!         two_body_2 = two_body_transcorr_factor(G1(2)%k, G1(4)%k)
! 
!         three_body_1 = three_body_transcorr_fac(hilbert_nI(:,1), G1(1)%k,G1(2)%k,G1(7)%k,-1)! + & 
!         three_body_2 = three_body_transcorr_fac(hilbert_nI(:,1),G1(2)%k,G1(1)%k,G1(4)%k,1)
! 
!         print *, "two-body: ", two_body_1, two_body_2
!         print *, "three-body: ", three_body_1, three_body_2
! 
!         print *, "excitation: (2,5) -> (7,8)" 
!         ! ex: should be + 
!         ! 2 5 
!         ! 7 8 
!         ! total 1 2 3 5 should be - 
! 
!         two_body_1 = two_body_transcorr_factor(G1(2)%k, G1(8)%k)! + & 
!         two_body_2 = two_body_transcorr_factor(G1(5)%k, G1(7)%k)
! 
!         three_body_1 = three_body_transcorr_fac(hilbert_nI(:,1), G1(2)%k,G1(5)%k,G1(8)%k,1)! + & 
!         three_body_2 = three_body_transcorr_fac(hilbert_nI(:,1),G1(5)%k,G1(2)%k,G1(7)%k,-1)
! 
!         print *, "two-body: ", two_body_1, two_body_2
!         print *, "three-body: ", three_body_1, three_body_2
! 
! 
!         print *, "excitation: (2,3) -> (6,7)" 
!         ! ex: should be - 
!         ! 2 3 
!         ! 6 7 
!         ! total 1 2 3 5 should be + 
! 
!         two_body_1 = two_body_transcorr_factor(G1(2)%k, G1(6)%k)! + & 
!         two_body_2 = two_body_transcorr_factor(G1(3)%k, G1(7)%k)
! 
!         three_body_1 = three_body_transcorr_fac(hilbert_nI(:,1), G1(2)%k,G1(3)%k,G1(6)%k,1)! + & 
!         three_body_2 = three_body_transcorr_fac(hilbert_nI(:,1),G1(3)%k,G1(2)%k,G1(7)%k,-1)
! 
!         print *, "two-body: ", two_body_1, two_body_2
!         print *, "three-body: ", three_body_1, three_body_2
! 
! 
!         print *, "excitation: (4,5) -> (1,8)" 
!         ! ex: should be - 
!         ! 4 5 
!         ! 1 8 
!         ! total 3 4 5 7 should be + 
! 
!         two_body_1 = two_body_transcorr_factor(G1(4)%k, G1(8)%k)! + & 
!         two_body_2 = two_body_transcorr_factor(G1(5)%k, G1(1)%k)
! 
!         three_body_1 = three_body_transcorr_fac(hilbert_nI(:,2), G1(4)%k,G1(5)%k,G1(8)%k,1)! + & 
!         three_body_2 = three_body_transcorr_fac(hilbert_nI(:,2),G1(5)%k,G1(4)%k,G1(1)%k,-1)
! 
!         print *, "two-body: ", two_body_1, two_body_2
!         print *, "three-body: ", three_body_1, three_body_2
! 
! 
! 
!         print *, "excitation: (3,4) -> (1,6)" 
!         ! ex: should be - 
!         ! 3 4 
!         ! 1 6 
!         ! total 3 4 5 7 should be - 
! 
!         two_body_1 = two_body_transcorr_factor(G1(4)%k, G1(6)%k)! + & 
!         two_body_2 = two_body_transcorr_factor(G1(3)%k, G1(1)%k)
! 
!         three_body_1 = three_body_transcorr_fac(hilbert_nI(:,2), G1(4)%k,G1(3)%k,G1(6)%k,1)! + & 
!         three_body_2 = three_body_transcorr_fac(hilbert_nI(:,2),G1(3)%k,G1(4)%k,G1(1)%k,-1)
! 
!         print *, "two-body: ", two_body_1, two_body_2
!         print *, "three-body: ", three_body_1, three_body_2
! 
! 
!         print *, "============== Excitation:  ======================"
!         print *, "nI: ", hilbert_nI(:,3)
!         print *, "nJ: ", hilbert_nI(:,4)
! 
!         print *, "H_ij: ", hamil_trancorr(3,4), hamil_trancorr(4,3)
! 
!         ! the excitation is from (3,8) -> (5,6)
!         print *, "excitation: (3,8) -> (5,6)"
!         ! ex: should have - sign 
!         ! 3 8 
!         ! 5 6
!         ! sign in total: 1 3 7 8 should be -
!         two_body = two_body_transcorr_factor(G1(3)%k, G1(5)%k) + & 
!                 two_body_transcorr_factor(G1(8)%k, G1(6)%k)
! 
!         print *, "two-body: ", two_body
! 
!         three_body_1 = three_body_transcorr_fac(hilbert_nI(:,3), G1(3)%k,G1(8)%k,G1(5)%k,-1)! + & 
!         three_body_2 = three_body_transcorr_fac(hilbert_nI(:,3),G1(8)%k,G1(3)%k,G1(6)%k,1)
! 
!         print *, "three-body: ", three_body_1, three_body_2

! 
!         print *, "and for k=0"
!         hilbert_nI(:,1) = [1,3,4,5]
!         hilbert_nI(:,2) = [1,2,3,7]
!         hilbert_nI(:,3) = [3,5,6,7]
!         hilbert_nI(:,4) = [1,5,7,8]
! 
!         ! first create the non-transcorrelated Hamiltonian
!         t_trans_corr_2body = .false. 
!         print *, "un-transcorrelated Hamiltonian: "
!         do i = 1, 4 
!             do j = 1, 4 
!                 hamil(i,j) = get_helement_k_space_hub(hilbert_nI(:,i),hilbert_nI(:,j))
!             end do
!             print *, hamil(i,:)
!         end do
!         
!         ! use the lapack routines to solve these quickly.. 
!         tmp_hamil = hamil
!         call dgeev('N','N',4,tmp_hamil,4,ev_real,ev_cmpl,left_ev,1,right_ev,1,work,12,info)
!         print *, "eigen-values: ", ev_real
! 
!         print *, "transcorrelated Hamiltonian: "
!         t_trans_corr_2body = .true.
!         do i = 1, 4 
!             do j = 1, 4 
!                 hamil_trancorr(i,j) = get_helement_k_space_hub(hilbert_nI(:,i),hilbert_nI(:,j))
!             end do
!             print *, hamil_trancorr(i,:)
!         end do
! 
!         tmp_hamil = hamil_trancorr
!         call dgeev('N','N',4,tmp_hamil,4,ev_real,ev_cmpl,left_ev,1,right_ev,1,work,12,info)
!         print *, "eigen-values: ", ev_real
! 
! !         call stop_all("here","for now")

    end subroutine test_4e_ms1

    subroutine setup_all(ptr, J, U)
        class(lattice), intent(in) :: ptr
        real(dp), intent(in), optional :: J, U

        nBasisMax = 0
        nullify(G1)
        nullify(tmat2d)
        deallocate(kPointToBasisFn)

        if (trim(ptr%get_name()) == 'tilted') then 
            ttilt = .true. 
        end if

        nBasis = 2*ptr%get_nsites()

        call setup_nbasismax(ptr)
        call setup_arr_brr(lat) 
        call setup_g1(ptr)
        call init_tmat_kspace(ptr)
!         call setup_tmat_k_space(ptr)
        call setup_kPointToBasisFn(ptr)
!         call setup_k_space_hub_sym(ptr)

        omega = real(ptr%get_nsites(),dp)

        bhub = -1.0

        if (present(U)) then 
            uhub = U 
        else 
            uhub = 1.0
        end if

        umat = h_cast(real(uhub,dp)/real(omega,dp))

        if (present(J)) then 
            trans_corr_param_2body = J 
        else 
            trans_corr_param_2body = 0.1
        end if

        three_body_prefac = 2.0_dp * (cosh(trans_corr_param_2body) - 1.0_dp) / real(omega**2,dp)

        ! after setup everything should be fine again.. or?
        ttilt = .false. 

    end subroutine setup_all

    subroutine exact_transcorrelation(lat_ptr, nI, J, U, hilbert_space_opt) 
        class(lattice), pointer, intent(in) :: lat_ptr 
        integer, intent(in) :: nI(:) 
        real(dp), intent(in) :: J(:), U 
        integer, intent(in), optional :: hilbert_space_opt(:,:)
        character(*), parameter :: this_routine = "exact_transcorrelation" 

        integer :: i, iunit, n_states, ind, k
        HElement_t(dp), allocatable :: hamil(:,:), hamil_trans(:,:), hamil_neci(:,:), &
                                       hamil_next(:,:), hamil_neci_next(:,:)
        real(dp), allocatable :: e_values(:), e_values_neci(:), e_vec(:,:), gs_vec(:)
        real(dp), allocatable :: e_vec_trans(:,:), t_mat_next(:,:), e_vec_next(:,:)
        real(dp) :: gs_energy, gs_energy_orig, hf_coeff_onsite(size(J))
        real(dp) :: hf_coeff_next(size(J)), hf_coeff_orig
        real(dp), allocatable :: neci_eval(:)
        integer, allocatable :: hilbert_space(:,:)
        character(30) :: filename, J_str
        
        ! then create the hilbert space 
        ! although make this an option to input it.. because i only 
        ! need to do that once actually.. 
        if (present(hilbert_space_opt)) then 
            ! if hilbert space is provided everything else should also be 
            ! setup already.. 
            allocate(hilbert_space(nel, size(hilbert_space_opt,2)), source = hilbert_space_opt)
            n_states = size(hilbert_space_opt,2)
        else
            call setup_system(lat_ptr, nI, J(1), U, hilbert_space) 
            n_states = size(hilbert_space,2) 
        end if

        print *, "total number of states: ", n_states 
        print *, "creating original hamiltonian: "
        t_trans_corr_2body = .false.
        hamil = create_hamiltonian(hilbert_space)

        print *, "diagonalizing original hamiltonian: " 
        allocate(e_values(n_states));        e_values = 0.0_dp
        allocate(e_vec(n_states, n_states)); e_vec = 0.0_dp
        allocate(gs_vec(n_states));          gs_vec = 0.0_dp
        do i = 1, size(hamil,1)
            do k = 1, size(hamil,2)
                if (isnan(hamil(i,k)) .or. is_inf(hamil(i,k))) print *, i,k,hamil(i,k)
            end do
        end do
        call eig(hamil, e_values, e_vec) 

        print *, "k-space hamiltonian: "
        call print_matrix(hamil)
        print *, "diagonal elements: e(k) + U/V"
        do i = 1, n_states 
            print *, hamil(i,i), "|", &
                sum(tmat2d(hilbert_space(:,i),hilbert_space(:,i))) + uhub/omega
        end do
        print *, "basis: " 
        print *, "i, k(1), k(2), k1 + k2, map(k1+k2)"
        do i = 1, n_states 
            print *,  hilbert_space(:,i), "|", &
                lat%get_k_vec(get_spatial(hilbert_space(1,i))), "|", & 
                lat%get_k_vec(get_spatial(hilbert_space(2,i))), "|", &
                lat%get_k_vec(get_spatial(hilbert_space(1,i))) + & 
                lat%get_k_vec(get_spatial(hilbert_space(2,i)))
        end do

        ! find the ground-state
        ind = minloc(e_values,1) 
        gs_energy_orig = e_values(ind) 

        ! how do i need to access the vectors to get the energy? 
        ! eigenvectors are stored in the columns!!
        gs_vec = abs(e_vec(:,ind))

        call sort(gs_vec)

        ! and flip order.. 
        gs_vec = gs_vec(n_states:1:-1)

        ! and i think i want the sorted by maximum of the GS
        print *, "original ground-state energy: ", gs_energy_orig
        ! and write the ground-state-vector to a file 
        iunit = get_free_unit()
        open(iunit, file = 'gs_vec_orig') 
        do i = 1, n_states
            write(iunit, *) gs_vec(i)
        end do
        close(iunit) 

        allocate(e_vec_trans(n_states, size(J)))
        e_vec_trans = 0.0_dp 
        allocate(e_vec_next(n_states, size(J)))
        e_vec_next = 0.0_dp

        hf_coeff_onsite = 0.0_dp
        hf_coeff_next = 0.0_dp

        ! also test that for the nearest neighbor transcorrelation 
        t_mat_next = get_tmat_next(lat_ptr, hilbert_space) 

        do i = 1, size(J) 
            print *, "J = ", J(i), ", U = ", U 

            write(J_str, *) J(i) 
            filename = 'gs_vec_trans_J_' // trim(adjustl((J_str)))

            trans_corr_param_2body = J(i)
            three_body_prefac = 2.0_dp * (cosh(trans_corr_param_2body) - 1.0_dp) / real(omega**2,dp)

            print *, "creating transformed hamiltonian: "
            hamil_trans = similarity_transform(hamil) 
            
            print *, "creating transformed hamiltonian with neighbor interaction" 
            hamil_next = similarity_transform(hamil, J(i) * t_mat_next)

            print *, "(and for testing purposes also create the neci-transcorrelated hamiltonian)"
            t_trans_corr_2body = .true.
            hamil_neci = create_hamiltonian(hilbert_space)
            t_trans_corr_2body = .false. 

            print *, "also test the the neighbor correlated neci hamiltonian" 
            
            t_trans_corr = .true. 
            trans_corr_param = J(i)*2.0
            hamil_neci_next = create_hamiltonian(hilbert_space) 
            t_trans_corr = .false. 

            neci_eval = calc_eigenvalues(hamil_neci)

            if (abs(gs_energy_orig - minval(neci_eval)) > 1.0e-12) then 
                print *, "original hamiltonian: "
                call print_matrix(hamil)
                print *, "on-site transcorr hamiltonian neci: "
                call print_matrix(hamil_neci)
                print *, "on-site transcorr transformed: "
                call print_matrix(hamil_trans)
                print *, "orig E0:    ", gs_energy_orig
                print *, "on-site E0: ", minval(neci_eval)

                call stop_all("here", "on-site transcorrelated energy not correct!")
            end if

            neci_eval = calc_eigenvalues(hamil_neci_next) 

            if (abs(gs_energy_orig - minval(neci_eval)) > 1.0e-12) then
                print *, "original hamiltonian: " 
                call print_matrix(hamil)
                print *, "neighbor transcorr neci: " 
                call print_matrix(hamil_neci_next) 
                print *, "neighbor transvorr transformed: " 
                call print_matrix(hamil_next) 
                print *, "orig E0:      ", gs_energy_orig
                print *, "next-site E0: ", minval(neci_eval)
                call stop_all("here", "neighbor transcorrelated energy not correct!")
            end if

            print *, "diagonalizing the transformed hamiltonian: " 
            call eig(hamil_trans, e_values, e_vec) 

            ! find the ground-state
            ind = minloc(e_values,1) 
            gs_energy = e_values(ind) 
            print *, "transformed ground-state energy: ", gs_energy 

            if (abs(gs_energy - gs_energy_orig) > 1.e-12) then 
                call stop_all("HERE!", "energy incorrect!")
            end if
            ! how do i need to access the vectors to get the energy? 
            gs_vec = abs(e_vec(:,ind))
            call sort(gs_vec)

            gs_vec = gs_vec(n_states:1:-1)

            hf_coeff_onsite(i) = gs_vec(1)

            e_vec_trans(:,i) = gs_vec
            ! and write the ground-state-vector to a file 
!             iunit = get_free_unit()
!             open(iunit, file = filename)
!             do k = 1, n_states
!                 write(iunit, *) gs_vec(k)
!             end do
!             close(iunit) 

            print *, "diagonalizing the neighbor transformed hamiltonian" 
            call eig(hamil_next, e_values, e_vec) 
            ind = minloc(e_values,1) 
            gs_energy = e_values(ind) 
            print *, "neighbor transformed ground-state energy: ", gs_energy 

            if (abs(gs_energy - gs_energy_orig) > 1.e-12) then 
                call stop_all("HERE!", "energy incorrect!")
            end if

            ! how do i need to access the vectors to get the energy? 
            gs_vec = abs(e_vec(:,ind))
            call sort(gs_vec)

            gs_vec = gs_vec(n_states:1:-1)

            hf_coeff_next(i) = gs_vec(1)

            e_vec_next(:,i) = gs_vec

         end do

        iunit = get_free_unit() 
        open(iunit, file = "hf_coeff_onsite")
        do i = 1, size(J)
            write(iunit, *) J(i), hf_coeff_onsite(i)
        end do

        iunit = get_free_unit()
        open(iunit, file = "hf_coeff_next")
        do i = 1, size(J) 
            write(iunit, *) J(i), hf_coeff_next(i)
        end do

        ! maybe plot all transformed into one file.. 
        iunit = get_free_unit() 
        open(iunit, file = "gs_vec_trans")
        ! the important quantitiy is J over U i guess or? 
        ! i am not sure.. 
        write(iunit, *) "# J: ", J
        do i = 1, n_states 
            write(iunit, *) e_vec_trans(i,:)
        end do
        close(iunit)

        open(iunit, file = "gs_vec_next")
        write(iunit, *) "# J: ", J 
        do i = 1, n_states
            write(iunit, *) e_vec_next(i,:) 
        end do

    end subroutine exact_transcorrelation

    function get_tmat_next(lat_ptr, hilbert_space) result(t_mat)
        ! in the k-space this essentially only is J*\sum_k \epsilon(k) n_k
        ! which is the setup tmat divided by 2
        class(lattice), pointer, intent(in) :: lat_ptr 
        integer, intent(in) :: hilbert_space(:,:) 
        real(dp) :: t_mat(size(hilbert_space,2),size(hilbert_space,2))

        integer :: i 

        t_mat = 0.0_dp 

        do i = 1, size(hilbert_space,2)
            t_mat(i,i) = sum(GetTMatEl(hilbert_space(:,i),hilbert_space(:,i))) / real(bhub,dp)
        end do



    end function get_tmat_next 

    subroutine test_general

        ! find the smallest system, where my code fails again.. 
        integer, allocatable :: nI(:), hilbert_nI(:,:) 
        integer(n_int), allocatable :: dummy(:,:)
        HElement_t(dp), allocatable :: hamil(:,:), hamil_trancorr(:,:), &
                                       trans_hamil(:,:), hamil_old(:,:)
        real(dp), allocatable :: eval(:), eval_neci(:), t_mat(:,:), evectors(:,:)
        integer :: n_states, iunit, n_pairs, i

        ! these are the quantitites to fix: 
        nOccAlpha = 4 
        nOccBeta = 4 
        nel = 8 
        lat => lattice('tilted', 2,2,1,.true.,.true.,.true.,'k-space')
        allocate(nI(nel)); nI = [1,2,3,4,5,6,9,10]

        call setup_all(lat)
        call setup_k_total(nI)

!         print *, "nBasisMax: "
!         do i = 1, size(nBasisMax,1)
!             print *, nBasisMax(i,:)
!         end do
!         print *, "G1: ", G1 
!         print *, "tmat: ", tmat2d

        ! i need a starting det
        call create_hilbert_space(nI, n_states, hilbert_nI, dummy, gen_all_excits_k_space_hubbard)

        print *, "n_states: ", n_states

        t_trans_corr_2body = .false. 
        hamil = create_hamiltonian(hilbert_nI)

        hamil_old = create_hamiltonian_old(hilbert_nI)

        allocate(eval(n_states))
        allocate(eval_neci(n_states))
        allocate(evectors(n_states,n_states))

!         eval = calc_eigenvalues(hamil)
        call eig(hamil, eval, evectors)

        call sort(eval)
        print *, "eigen-value orig: ", eval(1)


        t_trans_corr_2body = .true. 
        hamil_trancorr = create_hamiltonian(hilbert_nI)

        eval_neci = calc_eigenvalues(hamil_trancorr)
        call sort(eval_neci)
        print *, "eigen-value neci: ", eval_neci(1)

        print *, "diff: ", eval(1) - eval_neci(1) 

        allocate(t_mat(n_states,n_states))

        n_pairs = nOccAlpha * nOccBeta
        t_mat = get_tranformation_matrix(hamil, n_pairs) 

        trans_hamil = matmul(matmul(matrix_exponential(-t_mat),hamil),matrix_exponential(t_mat))

        eval = calc_eigenvalues(trans_hamil)

        call sort(eval) 
        print *, "eigen-value tran: ", eval(1)

        eval = calc_eigenvalues(hamil_old)
        call sort(eval) 
        print *, "eigen-value old:  ", eval(1)

        open(iunit,file='states')
        call print_matrix(transpose(hilbert_nI), iunit)
        close(iunit)

        open(iunit,file='hamil_orig')
        call print_matrix(hamil, iunit)
        close(iunit)

        open(iunit,file='hamil_neci')
        print *, "hamil neci: "
        call print_matrix(hamil_trancorr, iunit)
        close(iunit)

        print *, "==========================================="
        open(iunit,file='hamil_trans')
        print *, "hamil trans: "
        call print_matrix(trans_hamil, iunit)
        close(iunit)

        call stop_all("here", "now")

    end subroutine test_general

    subroutine test_8e_8orbs

        integer :: nI(8), n_states
        integer, allocatable :: hilbert_nI(:,:)
        integer(n_int), allocatable :: dummy(:,:)
        HElement_t(dp), allocatable :: hamil(:,:), hamil_trancorr(:,:), &
                                       trans_hamil(:,:)
        real(dp), allocatable :: eval(:), t_mat(:,:)

        nOccAlpha = 4
        nOccBeta = 4
        nel = 8 

        lat => lattice('tilted', 2, 2, 1, .true., .true., .true., 'k-space')

        call setup_all(lat)
        nI = [1,2,3,4,5,6,7,8]
        ! i need to set the momentum
        call setup_k_total(nI) 

        ! i need a starting det
        call create_hilbert_space(nI, n_states, hilbert_nI, dummy, gen_all_excits_k_space_hubbard)

        print *, "n_states: ", n_states

        t_trans_corr_2body = .false. 
        hamil = create_hamiltonian(hilbert_nI)

        allocate(eval(n_states))

        eval = calc_eigenvalues(hamil)

        call sort(eval)
        print *, "eigen-value orig: ", eval(1)


        t_trans_corr_2body = .true. 
        hamil_trancorr = create_hamiltonian(hilbert_nI)

        eval = calc_eigenvalues(hamil_trancorr)
        call sort(eval)
        print *, "eigen-value neci: ", eval(1)

        allocate(t_mat(n_states,n_states))

        t_mat = get_tranformation_matrix(hamil, 16) 

        trans_hamil = matmul(matmul(matrix_exponential(-t_mat),hamil),matrix_exponential(t_mat))

        eval = calc_eigenvalues(trans_hamil)

        call sort(eval) 
        print *, "eigen-value tran: ", eval(1)

        ! todo! ok, still a small mistake in the transcorrelated hamil!! 
        ! maybe thats why it behaves unexpected! 

!         call stop_all("here", "now")

    end subroutine test_8e_8orbs

    subroutine test_3e_ms1

        integer :: hilbert_nI(3,3), i, j, work(9), info
        HElement_t(dp) :: hamil(3,3), hamil_trancorr(3,3), tmp_hamil(3,3)
        real(dp) :: ev_real(3), ev_cmpl(3), left_ev(1,3), right_ev(1,3)
        real(dp) :: test(3,3), t_mat(3,3), trans_hamil(3,3)

        integer :: n_pairs

        nOccAlpha = 2
        nOccBeta = 1
        nel = 3

        ! i have to setup the whole system
        lat => lattice('chain',3,1,1,.true.,.true.,.true.,'k-space')

        call setup_all(lat)

        print *, "also test 3 electron system for consistency!" 

        hilbert_nI(:,1) = [2,3,4]
        hilbert_nI(:,2) = [4,5,6] 
        hilbert_nI(:,3) = [1,2,6]

        ! first create the non-transcorrelated Hamiltonian
        t_trans_corr_2body = .false. 
        print *, "un-transcorrelated Hamiltonian: "

        hamil = create_hamiltonian(hilbert_nI)
        call print_matrix(hamil)

        ! the originial hamiltonian also gives me the transformation matrix 
        ! do it the stupid way 
        t_mat = get_tranformation_matrix(hamil, 2) 
        
        ! use the lapack routines to solve these quickly.. 
        print *, "eigen-values: ", calc_eigenvalues(hamil)

        print *, "transcorrelated Hamiltonian: "
        t_trans_corr_2body = .true.

        hamil_trancorr = create_hamiltonian(hilbert_nI)
        call print_matrix(hamil_trancorr)

        print *, "eigen-values: ", calc_eigenvalues(hamil_trancorr)

        trans_hamil = matmul(matmul(matrix_exponential(-t_mat),hamil),matrix_exponential(t_mat))

        print *, "transformed hamiltonian: " 
        call print_matrix(trans_hamil)

        print *, "eigen-values: ", calc_eigenvalues(trans_hamil)

!         ! i know now, where there is an error-- between those 2 matrix 
!         ! elements: 
!         print *, "nI:", hilbert_nI(:,1)
!         print *, "nJ:", hilbert_nI(:,2)
!         print *, "H_ij: ", hamil_trancorr(1,2), hamil_trancorr(2,1)
! 
!         ! check the individual contribs here! 
!         ! excitation: (2,3) -> (5,6)
!         ! ex: should have a + sign
!         ! 2 3 
!         ! 5 6 
!         ! 2 3 4 -> should have an overall + sign
!         print *, "Excitation (2,3) -> (5,6)"
!         print *, "two-body: ", two_body_transcorr_factor(G1(2)%k,G1(6)%k), & 
!                                two_body_transcorr_factor(G1(3)%k,G1(5)%k)
! 
!         print *, "three_body: ", three_body_transcorr_fac(hilbert_nI(:,1), &
!                                     G1(2)%k,G1(3)%k,G1(6)%k,1), &
!                                  three_body_transcorr_fac(hilbert_nI(:,1), &
!                                     G1(3)%k,G1(2)%k,G1(5)%k,-1)
! 
!         print *, "excitation: (3,4) -> (1,6): "
!         ! ex: should have a - sign!
!         ! 3 4
!         ! 1 6
!         ! 2 3 4 -> should have an overal - sign
!         print *, "two-body: ", two_body_transcorr_factor(G1(3)%k,G1(1)%k), & 
!                                two_body_transcorr_factor(G1(4)%k,G1(6)%k)
! 
!         print *, "three_body: ", three_body_transcorr_fac(hilbert_nI(:,1), &
!                                     G1(3)%k,G1(4)%k,G1(1)%k,-1), &
!                                  three_body_transcorr_fac(hilbert_nI(:,1), &
!                                     G1(4)%k,G1(3)%k,G1(6)%k,1)
! 
!         print *, "excitation: (4,5) -> (1,2): "
!         ! ex: should have a + sign
!         ! 4 5
!         ! 1 2
!         ! 4 5 6 -> should have an overall + sign
!         print *, "two-body: ", two_body_transcorr_factor(G1(4)%k,G1(2)%k), & 
!                                two_body_transcorr_factor(G1(5)%k,G1(1)%k)
! 
!         print *, "three_body: ", three_body_transcorr_fac(hilbert_nI(:,2), &
!                                     G1(4)%k,G1(5)%k,G1(2)%k,1), &
!                                  three_body_transcorr_fac(hilbert_nI(:,2), &
!                                     G1(5)%k,G1(4)%k,G1(1)%k,-1)


        ! both sign conventions agree here! maybe it has to do with this?!


    end subroutine test_3e_ms1

    subroutine setup_g1_test

        class(lattice), pointer :: ptr 
        integer :: i 
        
        print *, "" 
        print *, "testing: setup_g1" 
        ptr => lattice('chain', 4, 1, 1,.true.,.true.,.true.,'k-space')
        nBasis = 8

        tpbc = .true. 
        treal = .false. 
        ttilt = .false.

        call setup_g1(ptr) 

        ! check ms: 
        do i = 1, 7, 2 
            call assert_equals(-1, G1(i)%ms)
        end do
        do i = 2, 8, 2
            call assert_equals(1, G1(i)%ms)
        end do

        call assert_equals([(0,i=1,8)], G1(1:8)%k(2),8)
        call assert_equals([(0,i=1,8)], G1(1:8)%k(3),8)
        call assert_equals([-1,-1,0,0,1,1,2,2],G1(1:8)%k(1),8)

        deallocate(G1)
        nBasisMax = 0
        nbasis = -1

    end subroutine setup_g1_test

    subroutine setup_nbasismax_test

        use SystemData, only: nBasisMax
        class(lattice), pointer :: ptr 
        print *, "" 
        print *, "testing: setup_nbasismax"
        ptr => lattice('chain', 4, 1, 1,.true.,.true.,.true.,'k-space')

        nBasis = 8
        tpbc = .true. 
        treal = .false. 

        call setup_nbasismax(ptr)

        call assert_equals(0, nBasisMax(1,3))
        call assert_equals(0, nBasisMax(1,4))
        call assert_equals(1, nBasisMax(2,4))
        call assert_equals(4, nBasisMax(1,5))
        call assert_equals(1, nBasisMax(2,5))
        call assert_equals(1, nBasisMax(3,5))

        nBasisMax = 0
        nbasis = -1

    end subroutine setup_nbasismax_test

    subroutine setup_tmat_k_space_test

        class(lattice), pointer :: ptr 
        integer :: i 

        ptr => lattice('chain', 4, 1, 1, .true., .true., .true.,'k-space')
        print *, "" 
        print *, "testing: setup_tmat_k_space"
        tpbc = .true. 
        treal = .false. 
        tOneElecDiag = .false. 
        ttilt = .false. 
        bhub = 1.0
        nBasis = 8 

        call setup_tmat_k_space(ptr)

        call assert_equals(h_cast(0.0_dp), GetTMatEl(1,2))

        call assert_equals(h_cast(2.0_dp*cos(-PI/2.0_dp)), GetTMatEl(1,1))
        call assert_equals(h_cast(2.0_dp*cos(-PI/2.0_dp)), GetTMatEl(2,2))
        call assert_equals(h_cast(2.0_dp*cos(PI/2.0_dp)), GetTMatEl(5,5))
        call assert_equals(h_cast(2.0_dp*cos(PI/2.0_dp)), GetTMatEl(6,6))
        call assert_equals(h_cast(2.0_dp), GetTMatEl(3,3))
        call assert_equals(h_cast(2.0_dp), GetTMatEl(4,4))
        call assert_equals(h_cast(-2.0_dp), GetTMatEl(7,7))
        call assert_equals(h_cast(-2.0_dp), GetTMatEl(8,8))

        deallocate(G1) 
        nBasisMax = 0
        deallocate(tmat2d)
        nbasis = -1

    end subroutine setup_tmat_k_space_test

    subroutine setup_kPointToBasisFn_test

        use SymExcitDataMod, only: kPointToBasisFn
        class(lattice), pointer :: ptr
        integer :: i
        print *, "" 
        print *, "testing: setup_kPointToBasisFn" 
        
        ptr => lattice('chain', 4, 1, 1, .true., .true., .true.,'k-space')

        call setup_kPointToBasisFn(ptr)

        call assert_equals(1, kPointToBasisFn(-1,0,0,1))
        call assert_equals(2, kPointToBasisFn(-1,0,0,2))
        call assert_equals(3, kPointToBasisFn(0,0,0,1))
        call assert_equals(4, kPointToBasisFn(0,0,0,2))
        call assert_equals(5, kPointToBasisFn(1,0,0,1))
        call assert_equals(6, kPointToBasisFn(1,0,0,2))
        call assert_equals(7, kPointToBasisFn(2,0,0,1))

        deallocate(kPointToBasisFn)

    end subroutine setup_kPointToBasisFn_test

    subroutine get_diag_helement_k_sp_hub_test

        print *, ""
        print *, "testing: get_diag_helement_k_sp_hub" 
        umat = 0.0_dp
        call assert_equals(h_cast(-4.0_dp), get_diag_helement_k_sp_hub([1,2,3,4]))

        umat = 2*uhub/omega
        call assert_equals(h_cast(-3.0_dp), get_diag_helement_k_sp_hub([1,2,3,4]))

        print *, "" 
        print *, "and now for 2-body transcorrelation: "
        t_trans_corr_2body = .true.
        ! test it for 0 transcorrelation
        trans_corr_param_2body = 0.0_dp
        call assert_equals(h_cast(-3.0_dp), get_diag_helement_k_sp_hub([1,2,3,4]))

        trans_corr_param_2body = 1.0_dp 

        three_body_prefac = 2.0_dp * (cosh(trans_corr_param_2body) - 1.0_dp) / real(omega**2,dp)

        umat = uhub/omega

        ! todo: tests for actual transcorrelation! 
        call assert_equals(h_cast(1.0_dp), get_diag_helement_k_sp_hub([1,2,3,4]))

    end subroutine get_diag_helement_k_sp_hub_test

    subroutine get_offdiag_helement_k_sp_hub_test

        integer :: ex(2,2)
        integer, allocatable :: nI(:)

        nel = 2 
        allocate(nI(nel)); nI = [3,4] 

        print *, ""
        print *, "testing: get_offdiag_helement_k_sp_hub" 
        t_trans_corr_2body = .false. 

        ! 0 due to spin-symmetry:
        ex(1,:) = [1,3]
        call assert_equals(h_cast(0.0_dp), get_offdiag_helement_k_sp_hub(nI,ex,.false.))
        ex(1,:) = [1,2]
        ex(2,:) = [6,4]
        call assert_equals(h_cast(0.0_dp), get_offdiag_helement_k_sp_hub(nI,ex,.false.))

        ! 0 due to momentum symmetry: 
        ex(2,:) = [3,4]
        call assert_equals(h_cast(0.0_dp), get_offdiag_helement_k_sp_hub(nI,ex,.false.))

        ! this should contribute:
        ex(1,:) = [1,6]
        ex(2,:) = [3,4]
        call assert_equals(h_cast(uhub/real(omega,dp)), get_offdiag_helement_k_sp_hub(nI,ex,.false.))
        ex(1,:) = [6,1]
        call assert_equals(h_cast(uhub/real(omega,dp)), get_offdiag_helement_k_sp_hub(nI,ex,.false.))

        ex(2,:) = [1,6]
        ex(1,:) = [3,4]
        call assert_equals(h_cast(-uhub/real(omega,dp)), get_offdiag_helement_k_sp_hub(nI,ex,.true.))

        t_trans_corr = .true. 
        trans_corr_param = 2.0_dp 

        call assert_equals(h_cast(uhub/real(omega,dp)*exp(4.0_dp)), get_offdiag_helement_k_sp_hub(nI,ex,.false.))

        ex(1,:) = [1,6]
        ex(2,:) = [3,4]

        call assert_equals(h_cast(uhub/real(omega,dp)*exp(-4.0_dp)), get_offdiag_helement_k_sp_hub(nI,ex,.false.))
        ex(1,:) = [6,1]
        call assert_equals(h_cast(uhub/real(omega,dp)*exp(-4.0_dp)), get_offdiag_helement_k_sp_hub(nI,ex,.false.))

        t_trans_corr = .false. 
        trans_corr_param = 0.0_dp 

        t_trans_corr_2body = .true. 
        trans_corr_param_2body = 1.0_dp

        ex(1,:) = [1,3]
        ex(2,:) = [2,4] 
        call assert_equals(h_cast(0.0_dp), get_offdiag_helement_k_sp_hub(nI,ex,.false.))

        ex(1,:) = [1,2]
        ex(2,:) = [3,4]
        call assert_equals(h_cast(0.0_dp), get_offdiag_helement_k_sp_hub(nI,ex,.false.))

        ! this should contribute:
        ex(1,:) = [1,6]
        ex(2,:) = [3,4]
        ! and have a triples contribution now!
        call assert_equals(h_cast(uhub/real(omega,dp)), get_offdiag_helement_k_sp_hub(nI,ex,.false.))
        ! this should be independent of order of electrons
        ex(1,:) = [6,1]
        call assert_equals(h_cast(uhub/real(omega,dp)), get_offdiag_helement_k_sp_hub(nI,ex,.false.))

!         three_body_prefac = 1.0_dp
        ex(1,:) = [1,6]
        ex(2,:) = [3,4]
        print *, "---------------------------------" 
        print *, "test order for matrix elements triples contrib to 'normal'"
        print *, "(1,6) -> (3,4): ", get_offdiag_helement_k_sp_hub(nI,ex,.false.)
        ex(1,:) = [6,1]
        print *, "(6,1) -> (3,4): ", get_offdiag_helement_k_sp_hub(nI,ex,.false.)
        ex(2,:) = [4,3]
        print *, "(6,1) -> (4,3): ", get_offdiag_helement_k_sp_hub(nI,ex,.false.)
        ex(1,:) = [1,6]
        print *, "(1,6) -> (4,3): ", get_offdiag_helement_k_sp_hub(nI,ex,.false.)
        print *, "---------------------------------" 

        ! 0 due to momentum symmetry
        ex(1,:) = [1,5]
        ex(2,:) = [3,7]
        call assert_equals(h_cast(0.0_dp), get_offdiag_helement_k_sp_hub(nI,ex,.false.))

        ex(1,:) = [2,6]
        ex(2,:) = [4,8]
        call assert_equals(h_cast(0.0_dp), get_offdiag_helement_k_sp_hub(nI,ex,.false.))

        ! this should contribute: 
        ex(1,:) = [1,3]
        ex(2,:) = [5,7] 
        call assert_equals(h_cast(-4.0_dp*three_body_prefac), get_offdiag_helement_k_sp_hub(nI,ex,.false.))

        ! the order in EX should not matter.. figure that out! 
        ex(1,:) = [3,1]
        call assert_equals(h_cast(-4.0_dp*three_body_prefac), get_offdiag_helement_k_sp_hub(nI,ex,.false.))

        print *," --------------------------" 
        print *, "test order of matrix elements for parallel excitations: " 
        ex(1,:) = [1,3]
        ex(2,:) = [5,7] 
        print *, "(1,3) -> (5,7): ", get_offdiag_helement_k_sp_hub(nI, ex,.false.)
        ex(1,:) = [3,1]
        print *, "(3,1) -> (5,7): ", get_offdiag_helement_k_sp_hub(nI, ex,.false.)
        ex(2,:) = [7,5]
        print *, "(3,1) -> (7,5): ", get_offdiag_helement_k_sp_hub(nI, ex,.false.)
        ex(1,:) = [1,3]
        print *, "(1,3) -> (7,5): ", get_offdiag_helement_k_sp_hub(nI, ex,.false.)
        print *," --------------------------" 

        
        ex(1,:) = [2,4]
        ex(2,:) = [6,8] 
        call assert_equals(h_cast(-4.0_dp*three_body_prefac), get_offdiag_helement_k_sp_hub(nI,ex,.false.))

        ! order of orbitals should also not matter! 
        ex(2,:) = [8,6] 
        call assert_equals(h_cast(-4.0_dp*three_body_prefac), get_offdiag_helement_k_sp_hub(nI,ex,.false.))
        ! and this should have opposite sign.... 
        ! why should this have the opposite sign?  the spin should not matter! 
        ex(2,:) = [1,3]
        ex(1,:) = [5,7] 
        call assert_equals(h_cast(-4.0_dp*three_body_prefac), get_offdiag_helement_k_sp_hub(nI,ex,.false.),1.0e-12)

        nel = 4 
        nI = [3,4,6,7]
        ex(1,:) = [4,6]
        ex(2,:) = [2,8]
        call assert_equals(h_cast(0.0_dp), get_offdiag_helement_k_sp_hub(nI, ex,.false.))

        ex(1,:) = [3,7] ! k = 2 
        ex(2,:) = [1,5] ! k = 0 
        call assert_equals(h_cast(0.0_dp), get_offdiag_helement_k_sp_hub(nI, ex,.false.))

        t_trans_corr_2body = .false. 

    end subroutine get_offdiag_helement_k_sp_hub_test

    subroutine get_helement_k_space_hub_test

        integer, allocatable :: ni(:), nJ(:)
        integer :: ex(2,3), ic_ret

        nel = 2
        allocate(nI(nel)); allocate(nJ(nel)); 

        print *, ""
        print *, "testing: get_helement_k_space_hub_test" 
        nI = [1,2]
        nJ = [3,4] 

        ic_ret = -1
        call assert_equals(h_cast(0.0_dp), get_helement_k_space_hub(nI,nJ,ic_ret))
        call assert_equals(2, ic_ret)

        call assert_equals(h_cast(uhub/omega), get_helement_k_space_hub([1,6],[3,4]))

        nel = 4
        deallocate(nI); deallocate(nJ); allocate(nJ(nel)); allocate(nI(nel)); 
        nI = [3,6,7,8]
        nJ = [1,2,5,8]
        ic_ret = -1
        call assert_equals(h_cast(0.0_dp), get_helement_k_space_hub(nI,nJ,ic_ret))
        call assert_equals(3, ic_ret) 

        t_trans_corr_2body = .true. 
        call assert_equals(h_cast(-4.0*three_body_prefac), get_helement_k_space_hub(nI,nJ))
        ! todo: more tests! 

        t_trans_corr_2body = .false. 

    end subroutine get_helement_k_space_hub_test

    subroutine pick_spin_opp_elecs_test

        integer, allocatable :: nI(:)
        integer :: elecs(2)
        real(dp) :: p_elec

        nel = 2 
        nOccBeta = 1
        nOccAlpha = 1
        allocate(nI(nel))

        print *, ""
        print *, "testing: pick_spin_opp_elecs"
        nI = [1,2] 
        call pick_spin_opp_elecs(nI, elecs, p_elec)

        call assert_equals(1.0_dp, p_elec) 
        if (elecs(1) == 1) then 
            call assert_equals(2, elecs(2)) 
        else if (elecs(1) == 2) then
            call assert_equals(1, elecs(2))
        end if

        nel = 4 
        nOccBeta = 2
        nOccAlpha = 2 
        deallocate(nI); allocate(nI(nel)); nI = [1,2,3,4]

        call pick_spin_opp_elecs(nI, elecs, p_elec) 
        call assert_equals(0.25_dp, p_elec) 
        call assert_true(.not. same_spin(nI(elecs(1)),nI(elecs(2))))

    end subroutine pick_spin_opp_elecs_test

    subroutine pick_from_cum_list_test

        integer :: ind
        real(dp) :: pgen
        print *, ""
        print *, "testing: pick_from_cum_list"
        call pick_from_cum_list([0.0_dp,1.0_dp],1.0_dp, ind, pgen)

        call assert_equals(2, ind) 
        call assert_equals(1.0_dp, pgen) 

        call pick_from_cum_list([1.0_dp,1.0_dp],1.0_dp, ind, pgen)

        call assert_equals(1, ind) 
        call assert_equals(1.0_dp, pgen) 

        call pick_from_cum_list([1.0_dp,2.0_dp],2.0_dp, ind, pgen)
        call assert_equals(0.5_dp, pgen) 


    end subroutine pick_from_cum_list_test

    subroutine create_ab_list_hubbard_test

        integer, allocatable :: nI(:), orb_list(:,:)
        integer(n_int), allocatable :: ilutI(:) 
        real(dp), allocatable :: cum_arr(:) 
        real(dp) :: cum_sum, cpt 
        integer :: tgt 

        nel = 4 
        allocate(nI(nel)); allocate(ilutI(0:niftot)); allocate(orb_list(8,2)); 
        allocate(cum_arr(8))

        print *, "" 
        print *, "testing: create_ab_list_hubbard"
        nI = [1,2,3,4]
        call EncodeBitDet(nI, ilutI)

        call create_ab_list_hubbard(nI, ilutI,[1,2], orb_list, cum_arr, cum_sum)

        call assert_equals(0.25_dp, cum_sum) 
        call create_ab_list_hubbard(nI, ilutI,[3,4], orb_list, cum_arr, cum_sum)
        call assert_equals(0.25_dp, cum_sum) 

        call create_ab_list_hubbard(nI, ilutI,[3,4], orb_list, cum_arr, cum_sum, 1, cpt)
        call assert_equals(0.25_dp, cum_sum) 
        call assert_equals(0.0_dp, cpt) 

        call create_ab_list_hubbard(nI, ilutI,[3,4], orb_list, cum_arr, cum_sum, 7, cpt)
        call assert_equals(0.5_dp, cpt)
        call create_ab_list_hubbard(nI, ilutI,[3,4], orb_list, cum_arr, cum_sum, 8, cpt)
        call assert_equals(0.5_dp, cpt)


        ! todo: i also have to do that for 2-body-transcorrelation, which 
        ! leads to parallel double excitations in the k-space hubbard 
        ! case -> check here if the get_orb_from_kpoints() function, works 
        ! correctly for ispn /= 2 and thub!

        ! todo: more tests! 

    end subroutine create_ab_list_hubbard_test

    subroutine calc_pgen_k_space_hubbard_test

        integer :: nI(4), ex(2,2)
        integer(n_int) :: ilutI(0:0)

        print *, "" 
        print *, "testing: calc_pgen_k_space_hubbard"

        ni = [1,2,3,4] 
        call EncodeBitDet(nI, ilutI) 

        ex(1,:) = [1,2]
        ex(2,:) = [5,6]

        call assert_equals(0.0_dp, calc_pgen_k_space_hubbard(nI, ilutI, ex, 1))
        call assert_equals(0.0_dp, calc_pgen_k_space_hubbard(nI, ilutI, ex, 0))
        call assert_equals(0.0_dp, calc_pgen_k_space_hubbard(nI, ilutI, ex, 3))
        call assert_equals(0.25_dp, calc_pgen_k_space_hubbard(nI, ilutI, ex, 2))

        ex(1,:) = [3,4]
        call assert_equals(0.0_dp, calc_pgen_k_space_hubbard(nI, ilutI, ex, 2))

        ex(2,:) = [7,8] 
        call assert_equals(0.25_dp, calc_pgen_k_space_hubbard(nI, ilutI, ex, 2))

    end subroutine calc_pgen_k_space_hubbard_test

    subroutine gen_excit_k_space_hub_test

        integer :: nI(4), ex(2,2), nJ(4)
        integer(n_int) :: ilutI(0:0), ilutJ(0:0)
        HElement_t(dp) :: hel
        real(dp) :: pgen 
        type(excit_gen_store_type) :: store
        logical :: tpar, t_found(6), found_all
        integer :: ic

        real(dp) :: p_elec = 0.25_dp

        nI = [1,2,3,4]
        call EncodeBitDet(nI, ilutI)

        print *, ""
        print *, "testing: gen_excit_k_space_hub" 

        t_found = .false. 
        found_all = .false. 

        do while (.not. found_all) 
            call gen_excit_k_space_hub(nI, ilutI, nJ, ilutJ, 0, ic, ex, & 
                    tpar, pgen, hel, store) 

            if (all(nJ == [3,4,5,6]) .and. .not. t_found(1)) then 
                t_found(1) = .true.
                ! do asserts: 
                call assert_equals(p_elec, pgen)
                call assert_true(.not. tpar)
            end if

            if (all(nJ == [2,3,5,8]) .and. .not. t_found(2)) then 
                t_found(2) = .true. 
                call assert_equals(p_elec/2.0, pgen)
                call assert_true(.not. tpar)
            end if

            if (all(nJ == [2,3,6,7]) .and. .not. t_found(3)) then 
                t_found(3) = .true. 
                call assert_equals(p_elec/2.0, pgen)
                call assert_true(.not. tpar)
            end if

            if (all(nJ == [1,4,5,8]) .and. .not. t_found(4)) then 
                t_found(4) = .true. 
                call assert_equals(p_elec/2.0, pgen)
                call assert_true(.not. tpar)
            end if

            if (all(nJ == [1,4,6,7]) .and. .not. t_found(5)) then 
                t_found(5) = .true. 
                call assert_equals(p_elec/2.0, pgen)
                call assert_true(.not. tpar)
            end if
            
            if (all(nJ == [1,2,7,8]) .and. .not. t_found(6))  then 
                t_found(6) = .true. 
                call assert_equals(p_elec, pgen)
                call assert_true(.not. tpar)
            end if

            found_all = all(t_found)

        end do

    end subroutine gen_excit_k_space_hub_test

    subroutine gen_parallel_double_hubbard_test

        integer :: nI(4), nJ(4), ex(2,2) 
        integer(n_int) :: ilutI(0:0), ilutJ(0:0) 
        real(dp) :: pgen 
        logical :: tpar, found_all, t_found(2)

        print *, "" 
        print *, "testing: gen_parallel_double_hubbard "

        t_trans_corr_2body = .true. 
        nI = [1,2,3,4]
        call EncodeBitDet(nI, ilutI)

        found_all = .false. 
        t_found = .false. 

        do while (.not. found_all)
            call gen_parallel_double_hubbard(nI, ilutI, nJ, ilutJ, ex, tPar, pgen) 

            if (all(nJ == [2,4,5,7]) .and. .not. t_found(1)) then 
                t_found(1) = .true. 
                call assert_equals(0.5_dp, pgen)
                call assert_true(tpar)
            end if

            if (all(nJ == [1,3,6,8]) .and. .not. t_found(2)) then 
                t_found(2) = .true. 
                call assert_equals(0.5_dp, pgen)
                call assert_true(tpar)
            end if 

            found_all = all(t_found)
        end do

        t_trans_corr_2body = .false. 

    end subroutine gen_parallel_double_hubbard_test

    subroutine gen_triple_hubbard_test

        integer :: nI(4), nJ(4), ex(2,3) 
        integer(n_int) :: ilutI(0:0), ilutJ(0:0) 
        real(dp) :: pgen 
        logical :: tpar, found_all, t_found(2)

        print *, "" 
        print *, "testing: gen_triple_hubbard "

        found_all = .false. 
        t_found = .false. 
        nI = [3,4,6,7]
        call EncodeBitDet(nI, ilutI) 

        t_trans_corr_2body = .true. 

        do while (.not. found_all) 
            call gen_triple_hubbard(nI, ilutI, nJ, ilutJ, ex, tpar, pgen) 

            if (all(nJ == [1,2,4,5]) .and. .not. t_found(1)) then 
                t_found(1) = .true. 
                call assert_equals(0.125_dp, pgen)
            end if
            found_all = t_found(1)

        end do

        nI = [3,4,5,8] 
        call EncodeBitDet(nI, ilutI) 

        found_all = .false. 
        t_found = .false. 

        do while (.not. found_all) 
            call gen_triple_hubbard(nI, ilutI, nJ, ilutJ, ex, tpar, pgen) 

            if (all(nJ == [1,2,3,6]) .and. .not. t_found(1)) then 
                t_found(1) = .true. 
                call assert_equals(0.125_dp, pgen)
            end if
            found_all = t_found(1)

        end do

        t_trans_corr_2body = .false.

    end subroutine gen_triple_hubbard_test

    subroutine pick_three_opp_elecs_test

        integer :: elecs(3), sum_ms
        real(dp) :: p_elec 

        nel = 3 
        nOccBeta = 2 
        nOccAlpha = 1 

        print *, ""
        print *, "testing: pick_three_opp_elecs"
        call pick_three_opp_elecs([1,2,3], elecs, p_elec, sum_ms)
        call assert_equals(1.0_dp, p_elec)
        call assert_equals(-1, sum_ms) 
        call assert_equals(6, sum(elecs))

        nOccAlpha = 2
        nOccBeta = 1 

        call pick_three_opp_elecs([1,2,4], elecs, p_elec, sum_ms)
        call assert_equals(1.0_dp, p_elec)
        call assert_equals(1, sum_ms) 
        call assert_equals(6, sum(elecs))

        nel = 5
        nOccAlpha = 4 
        call pick_three_opp_elecs([1,2,4,6,8], elecs, p_elec, sum_ms) 
        call assert_equals(1.0_dp/6.0_dp, p_elec) 
        call assert_equals(1, sum_ms) 
        call assert_true(any(elecs == 1))

        nOccBeta = 4 
        nOccAlpha = 1
        call pick_three_opp_elecs([1,3,5,7,8], elecs, p_elec, sum_ms) 
        call assert_equals(1.0_dp/6.0_dp, p_elec) 
        call assert_equals(-1, sum_ms) 
        call assert_true(any(elecs == 5))

        nel = 5
        nOccBeta = 3
        nOccAlpha = 2 
        call pick_three_opp_elecs([1,2,3,4,5], elecs, p_elec, sum_ms) 
        if (sum_ms == 1) then 
            call assert_equals(1.0_dp/10.0_dp, p_elec)
        else 
            call assert_equals(7.0_dp/60.0_dp, p_elec,1.0e-12)
        end if

        call pick_three_opp_elecs([1,2,3,4,5], elecs, p_elec, sum_ms) 
        if (sum_ms == 1) then 
            call assert_equals(1.0_dp/10.0_dp, p_elec)
        else 
            call assert_equals(7.0_dp/60.0_dp, p_elec, 1.0e-12)
        end if

        nel = 4
        nOccBeta = 2
        nOccAlpha = 2 
        call pick_three_opp_elecs([1,2,3,4], elecs, p_elec) 
        call assert_equals(0.25_dp, p_elec)

    end subroutine pick_three_opp_elecs_test

    subroutine pick_spin_par_elecs_test

        integer :: elecs(2), ispn
        real(dp) :: p_elec
        integer :: nI(6)
        
        print *, ""
        print *, "testing: pick_spin_par_elecs"
        nel = 2 
        nOccBeta = 2
        nOccAlpha = 0
        call pick_spin_par_elecs([1,3],elecs,p_elec, ispn)
        call assert_equals(1.0_dp, p_elec) 
        call assert_equals(1, ispn) 
        call assert_equals(3, sum(elecs))

        nOccAlpha = 2 
        nOccBeta = 0
        call pick_spin_par_elecs([2,4],elecs,p_elec, ispn)
        call assert_equals(1.0_dp, p_elec) 
        call assert_equals(3, ispn) 
        call assert_equals(3, sum(elecs))

        nel = 4 
        nOccBeta = 2
        call pick_spin_par_elecs([1,2,3,4], elecs, p_elec, ispn) 
        call assert_equals(0.5_dp, p_elec) 
        if (ispn == 1) then 
            call assert_equals(4, sum(elecs))
        else if (ispn == 3) then 
            call assert_equals(6, sum(elecs))
        end if
        
        nel = 6 
        nOccBeta = 3 
        nOccAlpha = 3

        call pick_spin_par_elecs([1,2,3,4,5,6], elecs, p_elec) 
        call assert_equals(1.0_dp/6.0_dp, p_elec) 
        nI = [1,2,3,4,5,6]
        call assert_true(same_spin(nI(elecs(1)),nI(elecs(2))))

        nel = 4
        nOccBeta = 2
        nOccAlpha = 2
    end subroutine pick_spin_par_elecs_test

    subroutine pick_a_orbital_hubbard_test

        integer(n_int) :: ilutI(0:0)
        integer :: orb, sum_ms
        real(dp) :: p_orb

        print *, ""
        print *, "testing: pick_a_orbital_hubbard "
        call EncodeBitDet([1,2,3,4], ilutI)

        call pick_a_orbital_hubbard(ilutI, orb, p_orb, -1) 
        call assert_true(orb == 6 .or. orb == 8) 
        call assert_equals(0.5_dp, p_orb)

        call pick_a_orbital_hubbard(ilutI, orb, p_orb, 1) 
        call assert_true(orb == 5 .or. orb == 7) 
        call assert_equals(0.5_dp, p_orb)

        call pick_a_orbital_hubbard(ilutI, orb, p_orb) 
        call assert_equals(0.25_dp, p_orb)

    end subroutine pick_a_orbital_hubbard_test

    subroutine pick_ab_orbitals_hubbard_test

        integer, allocatable :: nI(:)
        integer(n_int), allocatable :: ilutI(:) 
        integer :: orbs(2)
        real(dp) :: p_orb 

        allocate(nI(nel)); allocate(ilutI(0:niftot))

        print *, "" 
        print *, "testing: pick_ab_orbitals_hubbard"
        nI = [1,2,3,4]
        call EncodeBitDet(nI, ilutI)

        call pick_ab_orbitals_hubbard(nI, ilutI, [1,2], orbs, p_orb)

        call assert_equals(1.0_dp, p_orb) 
        call assert_equals(11, sum(orbs))

        call pick_ab_orbitals_hubbard(nI, ilutI, [3,4], orbs, p_orb)
        call assert_equals(1.0_dp, p_orb) 
        call assert_equals(15, sum(orbs))

    end subroutine pick_ab_orbitals_hubbard_test

    subroutine pick_bc_orbitals_hubbard_test

        integer:: nI(4)
        integer(n_int) :: ilutI(0:0) 
        integer :: orbs(2)
        real(dp) :: p_orb 


        t_trans_corr_2body = .true. 

        nI = [3,4,6,7]
        call EncodeBitDet(nI, ilutI)

        print *, ""
        print *, "testing: pick_bc_orbitals_hubbard"
        call pick_bc_orbitals_hubbard(nI, ilutI,[3,6,7],2,orbs,p_orb)
        call assert_equals(6, sum(orbs))
        call assert_equals(1.0_dp, p_orb) 

        nI = [3,4,5,8]
        call EncodeBitDet(nI, ilutI)
        call pick_bc_orbitals_hubbard(nI, ilutI,[4,5,8],1,orbs,p_orb)
        call assert_equals(8, sum(orbs))
        call assert_equals(1.0_dp, p_orb) 

    end subroutine pick_bc_orbitals_hubbard_test

    subroutine create_ab_list_par_hubbard_test

        integer:: nI(4), orb_list(4,2), tgt
        integer(n_int) :: ilutI(0:0) 
        real(dp) :: cum_sum, cpt, cum_arr(4)

        nI = [1,2,3,4]
        call EncodeBitDet(nI, ilutI) 

        t_trans_corr_2body = .true.
        print *, ""
        print *, "testing: create_ab_list_par_hubbard"
        call create_ab_list_par_hubbard(nI, ilutI, [1,3], orb_list, cum_arr, cum_sum)
        call assert_true(cum_sum > 0.0_dp)

        call create_ab_list_par_hubbard(nI, ilutI, [1,3], orb_list, cum_arr, cum_sum, 5, cpt)
        call assert_equals(0.5_dp, cpt) 
        call create_ab_list_par_hubbard(nI, ilutI, [1,3], orb_list, cum_arr, cum_sum, 7, cpt)
        call assert_equals(0.5_dp, cpt) 

        call create_ab_list_par_hubbard(nI, ilutI, [2,4], orb_list, cum_arr, cum_sum)
        call assert_true(cum_sum > 0.0_dp)

        call create_ab_list_par_hubbard(nI, ilutI, [2,4], orb_list, cum_arr, cum_sum, 1, cpt)
        call assert_equals(0.0_dp, cpt)
        call create_ab_list_par_hubbard(nI, ilutI, [2,4], orb_list, cum_arr, cum_sum, 5, cpt)
        call assert_equals(0.0_dp, cpt)

        call create_ab_list_par_hubbard(nI, ilutI, [2,4], orb_list, cum_arr, cum_sum, 6, cpt)
        call assert_equals(0.5_dp, cpt)

        nI = [1,2,4,5] 
        call EncodeBitDet(nI, ilutI) 
        call create_ab_list_par_hubbard(nI, ilutI, [1,5], orb_list, cum_arr, cum_sum)
        call assert_equals(0.0_dp, cum_sum)

        nI = [3,4,6,7]
        call EncodeBitDet(nI, ilutI)
        call create_ab_list_par_hubbard(nI, ilutI, [4,6], orb_list, cum_arr, cum_sum)

        print *, "" 
        print *, "cum_sum: ", cum_sum
        print *, "cum_arr: ", cum_arr
        print *, "orb_list(:,1):", orb_list(:,1)
        print *, "orb_list(:,2):", orb_list(:,2)

        t_trans_corr_2body = .false.

    end subroutine create_ab_list_par_hubbard_test

    subroutine pick_ab_orbitals_par_hubbard_test

        integer:: nI(4), orbs(2)
        integer(n_int) :: ilutI(0:0) 
        real(dp) :: p_orb 

        nI = [1,2,3,4]
        call EncodeBitDet(nI, ilutI) 

        t_trans_corr_2body = .true.

        print *, ""
        print *, "testing: pick_ab_orbitals_par_hubbard"
        call pick_ab_orbitals_par_hubbard(nI, ilutI, [1,3], orbs, p_orb)
        call assert_equals(12, sum(orbs))
        call assert_equals(1.0_dp, p_orb)

        call pick_ab_orbitals_par_hubbard(nI, ilutI, [2,4], orbs, p_orb)
        call assert_equals(14, sum(orbs))
        call assert_equals(1.0_dp, p_orb)

        nI = [1,2,4,5] 
        call EncodeBitDet(nI, ilutI) 
        call pick_ab_orbitals_par_hubbard(nI, ilutI, [1,5], orbs, p_orb)
        call assert_equals(0.0_dp, p_orb)

        t_trans_corr_2body = .false.

    end subroutine pick_ab_orbitals_par_hubbard_test

    subroutine get_transferred_momenta_test

        integer :: ex(2,2), k_vec_a(3), k_vec_b(3)

        print *, ""
        print *, "testing: get_transferred_momenta"
        ex(1,:) = [1,2]
        ex(2,:) = [3,4] 
        call get_transferred_momenta(ex, k_vec_a, k_vec_b) 

        call assert_equals(1, k_vec_a(1)) 
        call assert_equals(-1, k_vec_b(1))

        ex(1,:) = [1,3] 
        ex(2,:) = [5,7]
        call get_transferred_momenta(ex, k_vec_a, k_vec_b) 

        call assert_equals(2, k_vec_a(1))
        call assert_equals(-1, k_vec_b(1))

    end subroutine get_transferred_momenta_test

    subroutine create_bc_list_hubbard_test

        integer :: nI(4), orb_list(4,2), tgt 
        integer(n_int) :: ilutI(0:0) 
        real(dp) :: cum_arr(4), cum_sum, cpt 

        t_trans_corr_2body = .true. 
        print *, ""
        print *, "testing: create_bc_list_hubbard"
        nI = [1,2,3,4]
        call EncodeBitDet(nI, ilutI)
        call create_bc_list_hubbard(nI, ilutI, [1,2,3],6,orb_list, cum_arr, cum_sum)
        call assert_equals(0.0_dp, cum_sum) 

        call create_bc_list_hubbard(nI, ilutI, [1,2,4],5,orb_list, cum_arr, cum_sum)
        call assert_equals(0.0_dp, cum_sum) 

        nI = [3,4,6,7]
        call EncodeBitDet(nI, ilutI)
        call create_bc_list_hubbard(nI, ilutI, [3,6,7],2,orb_list, cum_arr, cum_sum)

        call assert_true(cum_sum > 0.0_dp) 

        nI = [3,4,5,8]
        call EncodeBitDet(nI, ilutI)
        call create_bc_list_hubbard(nI, ilutI, [4,5,8],1,orb_list, cum_arr, cum_sum)
        call assert_true(cum_sum > 0.0_dp) 
        call create_bc_list_hubbard(nI, ilutI, [4,5,8],1,orb_list, cum_arr, cum_sum, 4, cpt)
        call assert_equals(0.0_dp, cpt) 
        call create_bc_list_hubbard(nI, ilutI, [4,5,8],1,orb_list, cum_arr, cum_sum, 2, cpt)
        call assert_equals(0.5_dp, cpt)

        t_trans_corr_2body = .false. 

    end subroutine create_bc_list_hubbard_test

    subroutine get_3_body_helement_ks_hub_test

        integer :: nel, ex(2,3)
        integer, allocatable :: nI(:) 
        logical :: tpar

        tpar = .false.

        nel = 4 
        allocate(ni(nel))

        ni = [1,2,3,4]

        print *, "" 
        print *, "testing: get_3_body_helement_ks_hub" 
        ex(1,:) = [1,3,5] 
        ex(2,:) = [2,4,6]

        ! here spin does not fit:
        call assert_equals(h_cast(0.0_dp), get_3_body_helement_ks_hub(ni,ex,tpar))
        ex(1,:) = [1,2,3]
        ex(2,:) = [4,5,6]
        call assert_equals(h_cast(0.0_dp), get_3_body_helement_ks_hub(ni,ex,tpar))
        ex(2,:) = [5,6,7]

        ! and here momentum does not fit
        call assert_equals(h_cast(0.0_dp), get_3_body_helement_ks_hub(ni,ex,tpar))

        ! and here it should fit.
        ex(1,:) = [3,6,7]
        ex(2,:) = [1,2,5]
        call assert_equals(h_cast(-4*three_body_prefac), get_3_body_helement_ks_hub(ni,ex,tpar))

        ! and the order of the involved electrons should not change the 
        ! matrix element! ... damn.. it does.. i need to have some 
        ! convention i think.. as it is in the spin-opposite excitations 
        ! in the "normal" method.. 
        print *, "--------------------------" 
        print *, "testing order influence on sign: "
        print *, "(",ex(1,:),") -> (",ex(2,:),"): ", get_3_body_helement_ks_hub(nI, ex, .false.)
        ex(1,:) = [3,7,6]
        print *, "(",ex(1,:),") -> (",ex(2,:),"): ", get_3_body_helement_ks_hub(nI, ex, .false.)
        ex(1,:) = [6,3,7]
        print *, "(",ex(1,:),") -> (",ex(2,:),"): ", get_3_body_helement_ks_hub(nI, ex, .false.)
        ex(1,:) = [6,7,3] 
        print *, "(",ex(1,:),") -> (",ex(2,:),"): ", get_3_body_helement_ks_hub(nI, ex, .false.)
        ex(1,:) = [7,6,3]
        print *, "(",ex(1,:),") -> (",ex(2,:),"): ", get_3_body_helement_ks_hub(nI, ex, .false.)
        ex(2,:) = [1,5,2]
        print *, "(",ex(1,:),") -> (",ex(2,:),"): ", get_3_body_helement_ks_hub(nI, ex, .false.)
        ex(2,:) = [2,1,5]
        print *, "(",ex(1,:),") -> (",ex(2,:),"): ", get_3_body_helement_ks_hub(nI, ex, .false.)
        ex(2,:) = [5,1,2]
        print *, "(",ex(1,:),") -> (",ex(2,:),"): ", get_3_body_helement_ks_hub(nI, ex, .false.)
        ex(2,:) = [5,2,1]
        print *, "(",ex(1,:),") -> (",ex(2,:),"): ", get_3_body_helement_ks_hub(nI, ex, .false.)

        ! fix this sign incoherence above! 
        call assert_true(.false.)

    end subroutine get_3_body_helement_ks_hub_test

    subroutine check_momentum_sym_test

        print *, ""
        print *, "testing: check_momentum_sym"
        ! use the already setup up 4 site chain.. the input to this is 
        ! with spin-orbitals.. or no.. it is with spatial orbs! no it is 
        ! spin-orbital! but the spin is also checked for symmetry! 
        ! although it is not only momentum symmetry! it is also 
        ! spin symmetry!! 
        call assert_true(check_momentum_sym([1],[1]))
        call assert_true(.not.check_momentum_sym([1],[2]))
        call assert_true(check_momentum_sym([2],[2]))
        call assert_true(.not.check_momentum_sym([3],[1]))

        call assert_true(check_momentum_sym([1,2],[2,1]))
        call assert_true(check_momentum_sym([2,2],[2,2]))

        ! and it takes variable sizes of input..
        call assert_true(check_momentum_sym([1,5],[3,3]))

        call assert_true(check_momentum_sym([6,6],[4,8]))
        call assert_true(.not.check_momentum_sym([5,5],[4,8]))

        call assert_true(check_momentum_sym([5,8],[1,4]))
    end subroutine check_momentum_sym_test

    subroutine find_minority_spin_test

        print *, ""
        print *, "testing: find_minority_spin"
        call assert_equals(1, find_minority_spin([1,2,4]))
        call assert_equals(3, find_minority_spin([3,2,4]))

        call assert_equals(2, find_minority_spin([1,2,3]))
        call assert_equals(2, find_minority_spin([3,2,5]))

    end subroutine find_minority_spin_test

    subroutine calc_pgen_k_space_hubbard_transcorr_test

        integer :: nI(4), ex(2,3)
        integer(n_int) :: ilutI(0:0) 

        nI = [3,4,6,7]
        call EncodeBitDet(nI, ilutI)

        pDoubles = 0.8 
        pParallel = 0.2 

        t_trans_corr_2body = .true. 
        print *, "" 
        print *, "testing: calc_pgen_k_space_hubbard_transcorr"
        call assert_equals(0.0_dp, calc_pgen_k_space_hubbard_transcorr(ni,iluti,ex,0))
        call assert_equals(0.0_dp, calc_pgen_k_space_hubbard_transcorr(ni,iluti,ex,1))
        call assert_equals(0.0_dp, calc_pgen_k_space_hubbard_transcorr(ni,iluti,ex,4))
        
        ex(1,:) = [3,4,0] ! k = 0
        ex(2,:) = [2,5,0] ! k = 0
        ! this should contribute! 
        call assert_equals(0.8*0.8/4.0_dp, calc_pgen_k_space_hubbard_transcorr(nI,iluti,ex,2))

        ex(1,:) = [3,6,0] ! k = 1
        ex(2,:) = [1,8,0] ! k = 1
        call assert_equals(0.8*0.8/4.0_dp, calc_pgen_k_space_hubbard_transcorr(nI,iluti,ex,2))

        ex(1,:) = [3,7,0] ! k=2
        ex(2,:) = [1,5,0] ! k=0
        call assert_equals(0.0_dp, calc_pgen_k_space_hubbard_transcorr(ni,iluti,ex,2))

        ex(1,:) = [4,6,0] ! k = 1
        ex(2,:) = [2,8,0] ! k = 1
        ! should contribute! no, becuase diagonal part is 0!
        call assert_equals(0.0_dp, calc_pgen_k_space_hubbard_transcorr(ni,iluti,ex,2))

        nI = [1,4,6,7]
        call EncodeBitDet(nI, ilutI)
        call assert_equals(0.8*0.2/2.0_dp, calc_pgen_k_space_hubbard_transcorr(ni,iluti,ex,2))

        nI = [3,4,6,7]
        call EncodeBitDet(nI, ilutI)

        ! the triple should be: 
        ex(1,:) = [3,6,7] 
        ex(2,:) = [1,2,5] 
        call assert_equals(0.2/16.0_dp, calc_pgen_k_space_hubbard_transcorr(ni,iluti,ex,3),1.0e-12)

        t_trans_corr_2body = .false.

    end subroutine calc_pgen_k_space_hubbard_transcorr_test

    subroutine calc_pgen_k_space_hubbard_par_test

        integer :: nI(4), ex(2,2)
        integer(n_int) :: ilutI(0:0) 

        nI = [1,2,3,4] 
        call EncodeBitDet(nI, ilutI)

        t_trans_corr_2body = .true.
        print *, "" 
        print *, "testing: calc_pgen_k_space_hubbard_par"
        ex(1,:) = [1,3]
        ex(2,:) = [5,7] 

        call assert_equals(0.0_dp, calc_pgen_k_space_hubbard_par(nI,ilutI,ex,0))
        call assert_equals(0.0_dp, calc_pgen_k_space_hubbard_par(nI,ilutI,ex,1))
        call assert_equals(0.0_dp, calc_pgen_k_space_hubbard_par(nI,ilutI,ex,3))

        call assert_equals(0.5_dp, calc_pgen_k_space_hubbard_par(nI,ilutI,ex,2))
        ex(1,:) = [2,4]
        call assert_equals(0.0_dp, calc_pgen_k_space_hubbard_par(nI,ilutI,ex,2))
        ex(2,:) = [6,8]
        call assert_equals(0.5_dp, calc_pgen_k_space_hubbard_par(nI,ilutI,ex,2))
        ex(1,:) = [1,4]
        call assert_equals(0.0_dp, calc_pgen_k_space_hubbard_par(nI,ilutI,ex,2))

        t_trans_corr_2body = .false.

    end subroutine calc_pgen_k_space_hubbard_par_test

    subroutine calc_pgen_k_space_hubbard_triples_test

        integer :: nI(4), ex(2,3)
        integer(n_int) :: ilutI(0:0) 

        nI = [3,4,6,7]
        call EncodeBitDet(nI, ilutI) 

        t_trans_corr_2body = .true.
        print *, "" 
        print *, "testing: calc_pgen_k_space_hubbard_triples"
        ex(1,:) = [3,6,7]
        ex(2,:) = [1,2,5]

        call assert_equals(0.0_dp, calc_pgen_k_space_hubbard_triples(nI, ilutI, ex,0))
        call assert_equals(0.0_dp, calc_pgen_k_space_hubbard_triples(nI, ilutI, ex,1))
        call assert_equals(0.0_dp, calc_pgen_k_space_hubbard_triples(nI, ilutI, ex,2))

        call assert_equals(1.0_dp/(2*8.0_dp), calc_pgen_k_space_hubbard_triples(nI, ilutI, ex,3))

        ex(1,:) = [4,5,8] 
        call assert_equals(0.0_dp, calc_pgen_k_space_hubbard_triples(nI, ilutI, ex,3))
        ex(2,:) = [1,2,6]
        nI = [3,4,5,8] 
        call EncodeBitDet(nI, ilutI)
        call assert_equals(1.0_dp/(2*8.0_dp), calc_pgen_k_space_hubbard_triples(nI, ilutI, ex,3))

        t_trans_corr_2body = .false.

    end subroutine calc_pgen_k_space_hubbard_triples_test

    subroutine make_triple_test

        integer, allocatable :: nI(:), nJ(:)
        integer :: ex(2,3), ex2(2,3)
        logical :: tpar, tpar_2, tpar_3, tpar_4
        integer(n_int) :: ilutI(0:nifd), ilutJ(0:nifd)

        nel = 3 

        allocate(nI(nel))
        allocate(nJ(nel))

        print *, "" 
        print *, "testing: make_triple" 
        print *, "testing implicitly: FindExcitDet!"

        nI = [1,2,3] 
        call make_triple(nI,nJ,[1,2,3],[4,5,7],ex,tpar) 
        call assert_equals([4,5,7], nJ, 3)
        call assert_equals([1,2,3], ex(1,:),3)
        call assert_equals([4,5,7], ex(2,:),3)
        call assert_true(.not.tpar)

        ex2(1,1) = 3 
        call GetExcitation(nI,nJ,nel,ex2,tpar_4)
        call assert_equals(ex, ex2, 2,2)
        call assert_true(tpar .eqv. tpar_4)

        call EncodeBitDet(nI, ilutI)
        call EncodeBitDet(nJ, ilutJ)
        ex2(1,1) = 3 
        call GetBitExcitation(ilutI, ilutJ, ex2, tpar_3)
        call assert_equals(ex, ex2, 2, 3)

        call FindExcitDet(ex, nI, 3, tpar_2)
        call assert_true(tpar .eqv. tpar_2)
        call assert_true(tpar .eqv. tpar_3)

        ! and now more complicated stuff:
        nI = [1,2,3]
        call make_triple(nI,nJ,[1,3,2],[7,5,4],ex,tpar) 
        call assert_equals([4,5,7], nJ, 3)
        call assert_equals([1,2,3], ex(1,:),3)
        call assert_equals([4,5,7], ex(2,:),3)
        call assert_true(.not.tpar)

        ex2(1,1) = 3 
        call GetExcitation(nI,nJ,nel,ex2,tpar_4)
        call assert_equals(ex, ex2, 2,2)
        call assert_true(tpar .eqv. tpar_4)

        call EncodeBitDet(nJ, ilutJ)
        ex2(1,1) = 3 
        call GetBitExcitation(ilutI, ilutJ, ex2, tpar_3)
        call assert_equals(ex, ex2, 2, 3)

        call FindExcitDet(ex, nI, 3, tpar_2)
        call assert_true(tpar .eqv. tpar_2)
        call assert_true(tpar .eqv. tpar_3)

        nI = [1,2,5]
        call make_triple(nI,nJ,[3,1,2],[3,4,7],ex,tpar) 
        call assert_equals([3,4,7], nJ, 3)
        call assert_equals([1,2,5], ex(1,:),3)
        call assert_equals([3,4,7], ex(2,:),3)
        call assert_true(.not.tpar)

        ex2(1,1) = 3 
        call GetExcitation(nI,nJ,nel,ex2,tpar_4)
        call assert_equals(ex, ex2, 2,2)
        call assert_true(tpar .eqv. tpar_4)

        call EncodeBitDet(nI, ilutI)
        call EncodeBitDet(nJ, ilutJ)
        ex2(1,1) = 3 
        call GetBitExcitation(ilutI, ilutJ, ex2, tpar_3)
        call assert_equals(ex, ex2, 2, 3)

        ex(1,:) = [1,2,3]
        call FindExcitDet(ex, nI, 3, tpar_2)
        call assert_true(tpar .eqv. tpar_2)
        call assert_true(tpar .eqv. tpar_3)
        nI = [1,2,5]

        call make_triple(nI,nJ,[3,2,1],[8,7,3],ex,tpar) 
        call assert_equals([3,7,8], nJ, 3)
        call assert_equals([1,2,5], ex(1,:),3)
        call assert_equals([3,7,8], ex(2,:),3)
        call assert_true(.not.tpar)

        ex2(1,1) = 3 
        call GetExcitation(nI,nJ,nel,ex2,tpar_4)
        call assert_equals(ex, ex2, 2,2)
        call assert_true(tpar .eqv. tpar_4)

        call EncodeBitDet(nJ, ilutJ)
        ex2(1,1) = 3 
        call GetBitExcitation(ilutI, ilutJ, ex2, tpar_3)
        call assert_equals(ex, ex2, 2, 3)

        ex(1,:) = [1,2,3]
        call FindExcitDet(ex, nI, 3, tpar_2)
        call assert_true(tpar .eqv. tpar_2)
        call assert_true(tpar .eqv. tpar_3)
        nI = [1,2,5]

        call make_triple(nI,nJ,[3,2,1],[4,7,9],ex,tpar) 
        call assert_equals([4,7,9], nJ, 3)
        call assert_equals([1,2,5], ex(1,:),3)
        call assert_equals([4,7,9], ex(2,:),3)
        call assert_true(.not.tpar)
        ex(1,:) = [1,2,3]
        call FindExcitDet(ex, nI, 3, tpar_2)
        call assert_true(tpar .eqv. tpar_2)
        call assert_true(tpar .eqv. tpar_3)
        nI = [1,2,5]

        call EncodeBitDet(nJ, ilutJ)
        ex2(1,1) = 3 
        call GetBitExcitation(ilutI, ilutJ, ex2, tpar_3)
        call assert_equals(ex, ex2, 2, 3)

        call make_triple(nI, nJ, [1,2,3], [3,4,7], ex, tpar) 
        call assert_true(.not.tpar)

        ex2(1,1) = 3 
        call GetExcitation(nI,nJ,nel,ex2,tpar_4)
        call assert_equals(ex, ex2, 2,2)
        call assert_true(tpar .eqv. tpar_4)

        call EncodeBitDet(nJ, ilutJ)
        ex2(1,1) = 3 
        call GetBitExcitation(ilutI, ilutJ, ex2, tpar_3)
        call assert_equals(ex, ex2, 2, 3)

        ex(1,:) = [1,2,3]
        call FindExcitDet(ex, nI, 3, tpar_2)
        call assert_true(tpar .eqv. tpar_2)
        call assert_true(tpar .eqv. tpar_3)
        nI = [1,2,5]

        nel = 4 

        deallocate(nJ); allocate(NJ(nel))
        deallocate(nI); allocate(nI(nel))
        
        nI = [1,2,5,7]
        call make_triple(nI,nJ,[1,2,3],[3,6,9],ex,tpar)
        call assert_true(tpar)

        ex2(1,1) = 3 
        call GetExcitation(nI,nJ,nel,ex2,tpar_4)
        call assert_equals(ex, ex2, 2,2)
        call assert_true(tpar .eqv. tpar_4)

        call EncodeBitDet(nI, ilutI)
        call EncodeBitDet(nJ, ilutJ)
        ex2(1,1) = 3 
        call GetBitExcitation(ilutI, ilutJ, ex2, tpar_3)
        call assert_equals(ex, ex2, 2, 3)


        ex(1,:) = [1,2,3]
        call FindExcitDet(ex, nI, 3, tpar_2)
        call assert_true(tpar .eqv. tpar_2)
        call assert_true(tpar .eqv. tpar_3)

        nel = -1

    end subroutine make_triple_test

    subroutine make_double_test

        use get_excit, only: make_double
        integer, allocatable :: nJ(:),ni(:)
        integer :: ex(2,2), ex2(2,2)
        logical :: tpar, tpar_2, tpar_3, tpar_4
        integer(n_int) :: ilutI(0:nifd), ilutJ(0:nifd)

        print *, "" 
        print *, "testing: make_double" 
        print *, "to be consistent with the sign conventions! "

        ! to test this really strange sign convention also call all the other 
        ! routines here, which test sign.. 

        nel = 2
        allocate(nJ(nel))
        allocate(ni(nel)); 
        ni = [1,2]

        call make_double([1,2],nJ, 1,2, 3,4, ex, tpar)
        call assert_equals([3,4], nJ, 2) 
        call assert_equals([1,2],ex(1,:),2)
        call assert_equals(reshape([1,3,2,4],[2,2]),ex, 2,2)
        call assert_true(.not. tpar)

        call EncodeBitDet(nI, ilutI)
        call EncodeBitDet(nJ, ilutJ)
        ex2(1,1) = 2 
        call GetBitExcitation(ilutI, ilutJ, ex2, tpar_3)
        call assert_equals(ex, ex2, 2,2)

        ex2(1,1) = 2 
        call GetExcitation(nI,nJ,nel,ex2,tpar_4)
        call assert_equals(ex, ex2, 2,2)
        call assert_true(tpar .eqv. tpar_4)

        call FindExcitDet(ex, nI, 2, tpar_2)
        call assert_true(tpar .eqv. tpar_2)
        call assert_true(tpar .eqv. tpar_3)

        ni = [1,2]

        call make_double([1,2],nJ, 1,2, 5,4, ex, tpar)
        call assert_true(.not. tpar)
        call FindExcitDet(ex, nI, 2, tpar_2)
        call assert_true(tpar .eqv. tpar_2)
        ni = [1,2]

        ex2(1,1) = 2 
        call GetExcitation(nI,nJ,nel,ex2,tpar_4)
        call assert_equals(ex, ex2, 2,2)
        call assert_true(tpar .eqv. tpar_4)

        call EncodeBitDet(nJ, ilutJ)
        ex2(1,1) = 2 
        call GetBitExcitation(ilutI, ilutJ, ex2, tpar_3)
        call assert_equals(ex, ex2, 2,2)
        call assert_true(tpar .eqv. tpar_3)
        call assert_equals(ex, ex2, 2,2)

        call make_double([1,2],nJ, 1,2, 3,6, ex, tpar) 
        call assert_true(.not. tpar)
        call FindExcitDet(ex, nI, 2, tpar_2)
        call assert_true(tpar .eqv. tpar_2)
        ni = [1,2]

        ex2(1,1) = 2 
        call GetExcitation(nI,nJ,nel,ex2,tpar_4)
        call assert_equals(ex, ex2, 2,2)
        call assert_true(tpar .eqv. tpar_4)

        call EncodeBitDet(nJ, ilutJ)
        ex2(1,1) = 2 
        call GetBitExcitation(ilutI, ilutJ, ex2, tpar_3)
        call assert_equals(ex, ex2, 2,2)
        call assert_true(tpar .eqv. tpar_3)

        nel = 3 
        deallocate(nJ); allocate(nJ(nel))
        nI = [1,2,4]

        call make_double([1,2,4],nJ,1,2,5,6,ex,tpar)
        call assert_equals([4,5,6], nJ, 3)
        call assert_true(.not. tpar)

        ex2(1,1) = 2 
        call GetExcitation(nI,nJ,nel,ex2,tpar_4)
        call assert_equals(ex, ex2, 2,2)
        call assert_true(tpar .eqv. tpar_4)

        call FindExcitDet(ex, nI, 2, tpar_2)
        call assert_true(tpar .eqv. tpar_2)
        nI = [1,2,4]

        call EncodeBitDet(nI, ilutI)
        call EncodeBitDet(nJ, ilutJ)
        ex2(1,1) = 2 
        call GetBitExcitation(ilutI, ilutJ, ex2, tpar_3)
        call assert_equals(ex, ex2, 2,2)

        call make_double([1,2,4],nj,1,2,3,6,ex,tpar)
        call assert_true(tpar)

        ex2(1,1) = 2 
        call GetExcitation(nI,nJ,nel,ex2,tpar_4)
        call assert_equals(ex, ex2, 2,2)
        call assert_true(tpar .eqv. tpar_4)

        call FindExcitDet(ex, nI, 2, tpar_2)
        call assert_true(tpar .eqv. tpar_2)
        nI = [1,2,4]

        call EncodeBitDet(nJ, ilutJ)
        ex2(1,1) = 2 
        call GetBitExcitation(ilutI, ilutJ, ex2, tpar_3)
        call assert_equals(ex, ex2, 2,2)

        call make_double([1,2,4],nJ, 1, 2, 6, 7, ex, tpar)
        call assert_true(.not. tpar)
        nI = [1,2,4]

        ex2(1,1) = 2 
        call GetExcitation(nI,nJ,nel,ex2,tpar_4)
        call assert_equals(ex, ex2, 2,2)
        call assert_true(tpar .eqv. tpar_4)

        call FindExcitDet(ex, nI, 2, tpar_2)
        call assert_true(tpar .eqv. tpar_2)

        call EncodeBitDet(nJ, ilutJ)
        ex2(1,1) = 2 
        call GetBitExcitation(ilutI, ilutJ, ex2, tpar_3)
        call assert_equals(ex, ex2, 2,2)

        call make_double([1,2,3],nJ, 1, 2, 4, 7, ex, tpar)
        call assert_true(.not. tpar)
        nI = [1,2,3]

        ex2(1,1) = 2 
        call GetExcitation(nI,nJ,nel,ex2,tpar_4)
        call assert_equals(ex, ex2, 2,2)
        call assert_true(tpar .eqv. tpar_4)

        call EncodeBitDet(nI, ilutI)
        call EncodeBitDet(nJ, ilutJ)
        ex2(1,1) = 2 
        call GetBitExcitation(ilutI, ilutJ, ex2, tpar_3)
        call assert_equals(ex, ex2, 2,2)

        call FindExcitDet(ex, nI, 2, tpar_2)
        call assert_true(tpar .eqv. tpar_2)

        nel = -1


    end subroutine make_double_test

    subroutine three_body_transcorr_fac_test

        integer :: p(3), q(3), k(3)
        real(dp) :: test 

        nel = 4
        nOccBeta = 2
        nOccAlpha = 2

        print *, "" 
        print *, "testing: three_body_transcorr_fac"
        p = 0 
        q = 0
        k = 0

        call assert_equals(h_cast(0.0_dp), &
            three_body_transcorr_fac([1,2,3,4], p,q,k,1))
        call assert_equals(h_cast(0.0_dp), &
            three_body_transcorr_fac([1,2,3,4], p,q,k,-1))

        q(1) = 1 
        k(1) = 2 
        call assert_equals(h_cast(0.0_dp), &
            three_body_transcorr_fac([1,2,3,4], p,q,k,1),1e-12)

        p(1) = 2 
        q(1) = 0 
        k(1) = 1 
        call assert_equals(h_cast(0.0_dp), &
            three_body_transcorr_fac([1,2,3,4], p,q,k,1))

        ! is there a non-zero combination? 
        p(1) = 0
        q(1) = 2
        k(1) = 1
        call assert_equals(h_cast(three_body_prefac*8), & 
            three_body_transcorr_fac([1,2,3,4],p,q,k,1))
        call assert_equals(h_cast(three_body_prefac*8), & 
            three_body_transcorr_fac([1,2,3,4],p,q,k,-1))
        
        nOccBeta = 4 
        nOccAlpha = 0
        call assert_equals(h_cast(0.0_dp), & 
            three_body_transcorr_fac([1,3,5,7],p,q,k,-1))

        call assert_equals(h_cast(three_body_prefac*8), & 
            three_body_transcorr_fac([1,3,5,7],p,q,k,1))

        nOccBeta = 2
        nOccAlpha = 2

    end subroutine three_body_transcorr_fac_test

    subroutine two_body_transcorr_factor_test

        integer :: p(3), k(3)

        print *, "" 
        print *, "testing: two_body_transcorr_factor" 
        p = 0
        k = 0

        call assert_equals(h_cast(-4.0_dp/omega * (cosh(trans_corr_param_2body) - 1)), & 
            two_body_transcorr_factor(p,k))

        p(1) = 1
        call assert_equals(h_cast(0.0_dp), two_body_transcorr_factor(p,k),1.e-12)
        call assert_equals(h_cast(-4.0_dp/omega * sinh(trans_corr_param_2body)), & 
            two_body_transcorr_factor([2,0,0],[2,0,0]))

        call assert_equals(h_cast(4.0_dp/omega * sinh(trans_corr_param_2body)), & 
            two_body_transcorr_factor([0,0,0],[2,0,0]))

        call assert_equals(h_cast(-2.0_dp/omega *(exp(-trans_corr_param_2body)-1)), & 
            two_body_transcorr_factor([0,0,0],[1,0,0]),1.e-12)


        call assert_equals(h_cast(4.0_dp/omega * (cosh(trans_corr_param_2body) - 1)), & 
            two_body_transcorr_factor([2,0,0],[0,0,0]))

    end subroutine two_body_transcorr_factor_test

    subroutine epsilon_kvec_test

        integer :: k(3)

        ! depending on the lattice dimension.. 
        ! and also the tilted has other values or?? 
        print *, "" 
        print *, "testing: epsilon_kvec"
        k = 0 
        call assert_equals(h_cast(2.0_dp), epsilon_kvec(k))
        call assert_equals(h_cast(0.0_dp), epsilon_kvec([1,0,0]),1e-12)
        call assert_equals(h_cast(0.0_dp), epsilon_kvec([-1,0,0]),1e-12)
        call assert_equals(h_cast(-2.0_dp), epsilon_kvec([2,0,0]))

    end subroutine epsilon_kvec_test

    subroutine same_spin_transcorr_factor_test

        integer :: k(3) 
        integer, allocatable :: nI(:) 

        nel = 4
        allocate(nI(nel))

        nI = [1,2,3,4] 

        print *, "" 
        print *, "testing: same_spin_transcorr_factor" 

        call assert_equals(h_cast(three_body_prefac*4), same_spin_transcorr_factor([1,2,3,4],[0,0,0],1),1.e-12)
        call assert_equals(h_cast(three_body_prefac*4), same_spin_transcorr_factor([1,2,3,4],[0,0,0],-1),1.e-12)
        call assert_equals(h_cast(0.0_dp), same_spin_transcorr_factor([1,2,3,4],[1,0,0],-1),1.e-12)
        call assert_equals(h_cast(-three_body_prefac*4), same_spin_transcorr_factor([1,2,3,4],[2,0,0],-1),1.e-12)

        call assert_equals(h_cast(0.0_dp), same_spin_transcorr_factor([1,3,5,7],[0,0,0],1))
        call assert_equals(h_cast(0.0_dp), same_spin_transcorr_factor([1,3,5,7],[0,0,0],-1))

    end subroutine same_spin_transcorr_factor_test

    subroutine get_one_body_diag_test

        integer, allocatable :: nI(:)
        class(lattice), pointer :: ptr
        integer :: i 

        nel = 4

        allocate(nI(nel))
        nI = [1,2,3,4]
        
        ! i think i also want to use the lattice functionality in the 
        ! k-space hubbard model.. but i still have to think about how to do that!
        ! allow an additional input flag! 

!         ptr => lattice('chain', 4, 1, 1, .true., .true., .true.,'k-space')

        call stop_all("get_one_body_diag_test", "changed implementation!")
!         call setup_tmat_k_space(ptr) 
        print *, "" 
        print *, "testing: get_one_body_diag" 
        ! the spin = 1 means i want the diagonal contribution of the alpha 
        ! electrons! 
!         call assert_equals(h_cast(2.0_dp), get_one_body_diag(nI,1))
!         call assert_equals(h_cast(2.0_dp), get_one_body_diag(nI,-1))
!         call assert_equals(h_cast(4.0_dp), get_one_body_diag(nI))
! 
!         nI = [1,2,5,6]
!         call assert_equals(h_cast(0.0_dp), get_one_body_diag(nI,1),1.e-8)
!         call assert_equals(h_cast(0.0_dp), get_one_body_diag(nI,-1),1e-8)
!         call assert_equals(h_cast(0.0_dp), get_one_body_diag(nI),1.e-8)
! 
!         nI = [1,3,5,7]
!         call assert_equals(h_cast(0.0_dp), get_one_body_diag(nI,1))
!         call assert_equals(h_cast(0.0_dp), get_one_body_diag(nI,-1))
!         call assert_equals(h_cast(0.0_dp), get_one_body_diag(nI))


    end subroutine get_one_body_diag_test

    subroutine gen_excit_k_space_hub_test_stochastic

        integer :: nI(nel)

        print *, "" 
        print *, "try to implement a stochastic version to check the" 
        print *, "calculated pgen and the actual one.. " 
        ! i should write a general test-runner for this.. which takes 
        ! an excitation generator as an input.. that would be nice! 
        ! something like: " 
        ni = [1,2,3,4]
        
        call setup_k_total(nI) 

        call run_excit_gen_tester(gen_excit_k_space_hub, "gen_excit_k_space_hub", & 
            gen_all_excits=gen_all_excits_k_space_hubbard) 

    end subroutine gen_excit_k_space_hub_test_stochastic

    subroutine gen_excit_k_space_hub_transcorr_test_stoch

        integer :: nI(nel), n_excits, i, nJ(nel), n_triples
        integer(n_int), allocatable :: det_list(:,:)

        print *, "" 
        print *, "testing: gen_excit_k_space_hub_transcorr" 
        print *, "first for a system with no possible triples, due to momentum conservation"

        pDoubles = 0.8 
        pParallel = 0.2 
        t_trans_corr_2body = .true. 

        nI = [1,2,3,4] 
        call setup_k_total(nI)
        call gen_all_triples_k_space(nI, n_excits, det_list)

        print *, "number of triple excitations: ", n_excits
        do i = 1, n_excits
            call writebitdet(6, det_list(:,i),.true.)
        end do

        call gen_all_excits_k_space_hubbard(nI, n_excits, det_list)

        ! for this momentum sector there are no, triple excitations valid.. 
        ! so test that for now! 
        call run_excit_gen_tester(gen_excit_k_space_hub_transcorr_test, &
            "gen_excit_k_space_hub_transcorr",opt_ni = nI, & 
            gen_all_excits = gen_all_excits_k_space_hubbard) 

        do i = 1, n_excits 
            call decode_bit_det(nJ, det_list(:,i))
            call run_excit_gen_tester(gen_excit_k_space_hub_transcorr_test, &
                "gen_excit_k_space_hub_transcorr",opt_ni = nJ, & 
                gen_all_excits = gen_all_excits_k_space_hubbard) 
        end do


        print *, "" 
        print *, "and now for a system with triples: "
        nI = [3,4,6,7]
        call setup_k_total(nI) 

        call gen_all_triples_k_space(nI, n_triples, det_list)

        print *, "number of triple excitations: ", n_triples
        do i = 1, n_triples
            call writebitdet(6, det_list(:,i),.true.)
        end do

        call gen_all_excits_k_space_hubbard(nI, n_excits, det_list)

        ! for this momentum sector there are no, triple excitations valid.. 
        ! so test that for now! 
        call run_excit_gen_tester(gen_excit_k_space_hub_transcorr_test, &
            "gen_excit_k_space_hub_transcorr",opt_ni = nI, & 
            gen_all_excits = gen_all_excits_k_space_hubbard) 

        do i = 1, n_excits 
            call decode_bit_det(nJ, det_list(:,i))
            call run_excit_gen_tester(gen_excit_k_space_hub_transcorr_test, &
                "gen_excit_k_space_hub_transcorr",opt_ni = nJ, & 
                gen_all_excits = gen_all_excits_k_space_hubbard) 
        end do

    end subroutine gen_excit_k_space_hub_transcorr_test_stoch

end program test_k_space_hubbard
