program lowdin

    implicit none

    integer, parameter :: dp = selected_real_kind(15,307)
    integer :: i, j, k
    integer :: lwork, counter, naverage, nkeep
    integer :: ilen, jlen, nkeep_len
    integer :: npositive, info, ierr, nargs, stat
    integer :: hamil_unit, overlap_unit
    integer :: nrows, nrepeats
    real(dp), allocatable :: work(:)
    character(2) :: index_fmt, nkeep_fmt
    character(7) :: string_fmt
    character(25) :: ind1
    character(256) :: hamil_filename, overlap_filename
    character(256) :: nrows_str, nrepeats_str, temp_str
    character(len=*), parameter :: stem = "lowdin"
    logical :: does_exist

    real(dp), allocatable :: overlap_estimates(:)
    real(dp), allocatable :: hamil_estimates(:)

    real(dp), allocatable :: hamil_mean(:,:)
    real(dp), allocatable :: overlap_mean(:,:)
    real(dp), allocatable :: overlap_eigenvecs(:,:)
    real(dp), allocatable :: transform_matrix(:,:)
    real(dp), allocatable :: inter_hamil(:,:)
    real(dp), allocatable :: final_hamil(:,:)
    real(dp), allocatable :: eigenvec_krylov(:,:)
    real(dp), allocatable :: hamil_eigv(:)
    real(dp), allocatable :: init_overlaps(:)
    real(dp), allocatable :: overlap_eigv(:)

    nargs = command_argument_count()
    if (nargs < 4) then
        write(6,'()') 'Please provide exactly 4 input arguments. The first two &
                       &should be the number of rows in the Hamiltonian and &
                       &then the number of repeats to consider. The third &
                       &should be the name of the Hamiltonian file, and the &
                       &fourth the name of the overlap matrix file.'
        stop 999
    end if

    call get_command_argument(1,nrows_str)
    read(nrows_str, '(i10)') nrows

    call get_command_argument(2,nrepeats_str)
    read(nrepeats_str, '(i10)') nrepeats

    call get_command_argument(3,hamil_filename)
    inquire(file=hamil_filename, exist=does_exist)
    if (.not. does_exist) then
        write(6,'(a21)') 'hamil file '//trim(hamil_filename)//' does not exist.'
        stop 999
    end if

    call get_command_argument(4,overlap_filename)
    inquire(file=overlap_filename, exist=does_exist)
    if (.not. does_exist) then
        write(6,'(a21)') 'hamil file '//trim(overlap_filename)//' does not exist.'
        stop 999
    end if

    allocate(hamil_estimates(nrepeats))
    allocate(overlap_estimates(nrepeats))

    allocate(hamil_mean(nrows, nrows))
    allocate(overlap_mean(nrows, nrows))
    allocate(overlap_eigenvecs(nrows, nrows))
    allocate(overlap_eigv(nrows))

    hamil_estimates = 0.0_dp
    hamil_mean = 0.0_dp
    overlap_estimates = 0.0_dp
    overlap_mean = 0.0_dp

    hamil_unit = 10
    open(hamil_unit, file=hamil_filename, status='old')
    ! Read in the hamil elements and errors to the above arrays.
    do i = 1, nrows
        ilen = ceiling(log10(real(abs(i)+1)))
        do j = i, nrows
            jlen = ceiling(log10(real(abs(j)+1)))
            write(index_fmt,'(a1,i1)') "a", ilen + jlen + 3

            read(hamil_unit, '('//index_fmt//','//trim(nrepeats_str)//'(1x,es19.12))', iostat=stat) &
                temp_str, hamil_estimates

            do k = 1, nrepeats
                hamil_mean(i,j) = hamil_mean(i,j) + hamil_estimates(k)
            end do
            hamil_mean(i,j) = hamil_mean(i,j)/nrepeats
            hamil_mean(j,i) = hamil_mean(i,j)
        end do
    end do
    close(hamil_unit)

    overlap_unit = 10
    open(overlap_unit, file=overlap_filename, status='old')
    do i = 1, nrows
        ilen = ceiling(log10(real(abs(i)+1)))
        do j = i, nrows
            jlen = ceiling(log10(real(abs(j)+1)))
            write(index_fmt,'(a1,i1)') "a", ilen + jlen + 3

            read(hamil_unit, '('//index_fmt//','//trim(nrepeats_str)//'(1x,es19.12))', iostat=stat) &
                temp_str, overlap_estimates

            write(*,*) overlap_estimates

            do k = 1, nrepeats
                overlap_mean(i,j) = overlap_mean(i,j) + overlap_estimates(k)
            end do
            overlap_mean(i,j) = overlap_mean(i,j)/nrepeats
            overlap_mean(j,i) = overlap_mean(i,j)
        end do
    end do
    close(overlap_unit)

    ! Create the workspace for the diagonaliser.
    lwork = max(1,3*nrows-1)
    allocate(work(lwork), stat=ierr)

    overlap_eigenvecs = overlap_mean

    ! Now perform the diagonalisation.
    call dsyev('V', 'U', nrows, overlap_eigenvecs, nrows, overlap_eigv, work, lwork, info)

    npositive = 0
    write(*,'(4("-"),a26,40("-"))') "Overlap matrix eigenvalues"
    do i = 1, nrows
        write(*,'(1x,es19.12)') overlap_eigv(i)
        if (overlap_eigv(i) > 0.0_8) npositive = npositive + 1
    end do

    do nkeep = 1, npositive

        allocate(transform_matrix(1:nrows, 1:nkeep))
        allocate(inter_hamil(1:nrows, 1:nkeep))
        allocate(final_hamil(1:nkeep, 1:nkeep))
        allocate(eigenvec_krylov(1:nrows, 1:nkeep))
        allocate(hamil_eigv(1:nkeep))
        allocate(init_overlaps(1:nkeep))

        counter = 0
        do i = nrows-nkeep+1, nrows
            counter = counter + 1
            transform_matrix(:,counter) = overlap_eigenvecs(:, i)/sqrt(overlap_eigv(i))
        end do

        inter_hamil = matmul(hamil_mean, transform_matrix)
        final_hamil = matmul(transpose(transform_matrix), inter_hamil)

        call dsyev('V', 'U', nkeep, final_hamil, nkeep, hamil_eigv, work, lwork, info)

        eigenvec_krylov = matmul(transform_matrix, final_hamil)
        init_overlaps = matmul(overlap_mean(1,:), eigenvec_krylov)

        nkeep_len = ceiling(log10(real(abs(nkeep)+1)))
        write(nkeep_fmt,'(a1,i1)') "i", nkeep_len
        write(string_fmt,'(i2,a5)') 15-nkeep_len, '("-")'
        write(*,'(/,4("-"),a37,1x,'//nkeep_fmt//',1x,a12,'//string_fmt//')') &
            "Eigenvalues and overlaps when keeping", nkeep, "eigenvectors"
        do i = 1, nkeep
            write(*,'(1x,es19.12,1x,es19.12)') hamil_eigv(i), init_overlaps(i)
        end do

        deallocate(transform_matrix)
        deallocate(inter_hamil)
        deallocate(final_hamil)
        deallocate(eigenvec_krylov)
        deallocate(hamil_eigv)
        deallocate(init_overlaps)

    end do

    deallocate(work)

    deallocate(hamil_estimates)
    deallocate(overlap_estimates)

    deallocate(hamil_mean)
    deallocate(overlap_mean)
    deallocate(overlap_eigenvecs)

end program lowdin
