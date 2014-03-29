program main
! This program reads in an FCIDMP input file and performs rotation of the 
! orbitals in order to break symmetry and maximally localise or delocalise 
! the rotated orbitals.
! The rotations are performed pair-wise by calling the Rotate2Orbs routine.
! compile with: f95 -F2008 -o fcidump_rotation.x fcidump_rotation.f90 -llapack -lblas -latlas
! / -lacml

implicit none

! constant date

!integer, parameter :: sp = selected_real_kind(6,37)
integer, parameter :: dp = selected_real_kind(15,307)
!integer, parameter :: qp = selected_real_kind(33,4931)
!integer, parameter :: int32 = selected_int_kind(8)
integer, parameter :: int64 = selected_int_kind(15)

real(dp), parameter :: pi = 3.1415926535897932384626433832795028841971693993751_dp

integer :: norb,nelec,ms2,isym,syml(1000),iuhf,symlz(1000),nprop(3)
integer :: propbitlen
logical :: uhf
real(dp), allocatable :: umat(:,:,:,:),tmat(:,:),arr(:),temp_inds(:,:)
real(dp), allocatable :: temp_tmat(:,:)
real(dp) :: ecore,diff
complex(dp) :: cz
real(dp) :: z
integer(int64) :: orbsym(1000)
integer :: i,j,k,l,ierr,ispin,m,nsets
integer :: specifier,numb,loc,l1,l2,l3,l4,l5
real(dp), allocatable :: selfint(:)
real(dp) :: selfint_old,selfint_prev,selfint_curr
real(dp),allocatable :: transform(:,:),fdiag(:)
real(dp),allocatable :: trans_2orbs_coeffs(:,:)
integer :: startindex,iumat,jumat,rotatepairs,iter
integer :: indices(1000)
integer :: scheme
integer, allocatable :: energyorder(:)
logical :: exists,trealinds,tmolpro,trotdegen,complexint
logical, allocatable :: localdelocal(:)
integer, allocatable :: sets(:),rotate_list(:,:),pairlist(:,:)
namelist /fci/ norb,nelec,ms2,orbsym,isym,iuhf,uhf,syml,symlz,propbitlen,nprop
    
    ! Defaults for rot_params input file
    trealinds = .true.
    tmolpro = .false.
    trotdegen = .true.
    specifier = 1
    loc = 0
    startindex = 0
    numb = 2
    indices(:) = 0

    ! Remaining defaults 
    uhf = .false.
    propbitlen = 0
    nprop = 0
    iuhf = 0

    write(6,*) '---------------------------------------------------'
    write(6,*) 'Rotating Orbitals from an FCIDUMP file'
    write(6,*) '---------------------------------------------------'
 
    ! Read in input options for rotation
    inquire(file='rot_params',exist=exists)
    if (.not.exists) then
        stop 'rot_params input file does not exist'
    endif
    write(6,*) 'Reading in rot_params input file...'
    open(9,file='rot_params',status='old',action='read')
    read(9,*) scheme         ! scheme: specifies which rotation scheme should be used; if scheme=1 an Edminston-Ruedenberg type localisation is employed which maximises (minimises the self-interactions energy <ii|ii>, if scheme=2 the quantity 
                             !\sum_{k=1,2}\sum_{j in occ} <kj|kj>-<kj|jk> is maximised (minimised)
    read(9,*) trealinds      ! integrals in FCIDUMP file are real (if .true.) and complex 
                             ! (if .false.)
    read(9,*) tmolpro        ! Molpro FCIDUMP file (if .true.) or Qchem FCIDUMP file 
                             ! (if .false.)
    read(9,*) trotdegen      ! sets of degenerate orbitals will also be rotated amongst 
                             ! each other (if.true.) and otherwise not (if.false.)
    read(9,*) diff           ! criterion for self-consistency of rotations, diff is the difference between the selfinteractions of two subsequent rotation cycles, if the difference between two subsequent rotation cycles is smaller than diff, the rotation is deemed to be self-consistent and stopped
    read(9,*) specifier      ! specifier: indication how the sets of orbitals to be rotated are specified: if 1: the following is just a list of the indices of the orbitals to rotate; if 2: the following is a list of indices of the starting and beginning indices of the orbitals to rotate
    close(9)

    if (scheme.eq.1) then
        write(6,*) '------------------------------------------------------------------'
        write(6,*) 'Edminston-Ruedenberg type procedure using self-interaction <ii|ii>'
        write(6,*) '------------------------------------------------------------------'
    elseif (scheme.eq.2) then
        write(6,*) '-------------------------------------------------------------'
        write(6,*) 'Procedure uses \sum_{k=1,2} \sum_{j in occ} (<kj|kj>-<kj|jk>)'
        write(6,*) '-------------------------------------------------------------'
    endif

   ! Read in original FCIDUMP
    inquire(file='FCIDUMP',exist=exists)
    if (.not.exists) then
        stop 'FCIDUMP file does not exist'
    endif
    
    write(6,*) 'Reading in FCIDUMP file...'
    open(8,file='FCIDUMP',status='old',form='formatted',action='read')
    ! Read in all symmetry labels although most of this will be 
    ! discarded
    read(8,fci)

    ! Molpro's FCIDUMP files always indicate the number of spatial orbital
    ! The number of spin orbitals is thus
    if (uhf.and.tmolpro) then
        norb = 2*norb
    endif

    ! not using any symmetry for storing integrals
    allocate(umat(norb,norb,norb,norb),stat=ierr)
    if (ierr.ne.0) then
        stop 'Error allocating Umat'
    endif
    allocate(tmat(norb,norb),stat=ierr)
    if (ierr.ne.0) then
        stop 'Error allocating Tmat'
    endif
    allocate(arr(norb),stat=ierr)
    if (ierr.ne.0) then
        stop 'Error allocating Arr'
    endif
    allocate(localdelocal(norb),stat=ierr)
    if (ierr.ne.0) then
        stop 'Error allocating localisation vector'
    endif
    allocate(sets(norb),stat=ierr)
    if (ierr.ne.0) then
        stop 'Error allocating vector with sets'
    endif
    allocate(rotate_list(norb,norb),stat=ierr)
    if (ierr.ne.0) then
        stop 'Error allocating list of orbitals to rotate'
    endif
    allocate(selfint(norb),stat=ierr)
    if (ierr.ne.0) then
        stop 'Error allocating list for selfiteractions'
    endif
    allocate(transform(norb,norb),stat=ierr)
    if (ierr.ne.0) then
        stop 'Error allocating transformation matrix'
    endif
    allocate(trans_2orbs_coeffs(2,2),stat=ierr)
    if (ierr.ne.0) then
        stop 'Error allocating 2 orbital transformation matrix'
    endif
    allocate(temp_inds(norb,norb),stat=ierr)
    if (ierr.ne.0) then
        stop 'Error allocating temporary Umat'
    endif
    allocate(temp_tmat(norb,norb),stat=ierr)
    if (ierr.ne.0) then
        stop 'Error allocating temporary tmat'
    endif
    allocate(fdiag(norb),stat=ierr)
    if (ierr.ne.0) then
        stop 'Error allocating diagonal fock matrix elements'
    endif
    allocate(energyorder(norb),stat=ierr)
    if (ierr.ne.0) then
        stop 'Error allocating list for ordered orbitals'
    endif



    umat(:,:,:,:) = 0.0_dp
    tmat(:,:) = 0.0_dp
    arr(:) = 0.0_dp
    ecore = 0.0_dp
    localdelocal(:) = .true.
    sets(:) = 0
    rotate_list(:,:) = 0
    selfint(:) = 0.0_dp
    transform(:,:) = 0.0_dp
    temp_inds(:,:) = 0.0_dp

    ! Umat contains integrals in physical notation
    ! <ij|kli>
    ! FCIDUMP files are in chemical notation, i.e.
    ! (ik|jl) = <ij|kl>
    ! real integrals
    ! if UHF FCIDUMP files the orbitals are ordered 
    ! firstly by spatial indices i and then by spin a/b, i.e.
    ! (i,a),(i,b),(j,a),(j,b) ...
    if (trealinds) then
        complexint = .false.
        ispin = 1
        do
            if (.not.tmolpro) then 
                read(8,'(1X,G20.12,4I3)',end=199) z,i,k,j,l
            elseif (tmolpro) then
                ! Molpro writes out integrals to a greater precision
                read(8,*,end=199) z,i,k,j,l
            endif
            if (uhf.and.tmolpro) then
                ! Molpro writes out spatial orbitals indices which need 
                ! to be transfered to spin orbitals indices
                ! UHF FCIDUMP files in Molpro are written out as
                !1: aaaa
                !2: bbbb
                !3: aabb
                !4: aa
                !5: bb
                !with delimiters of 0.000000 0 0 0 0
                ! need to transform into spin orbital indices
                if ((i.eq.0).and.(k.eq.0).and.(j.eq.0).and.(l.eq.0).and.&
                    &(ispin.ne.6)) then
                    if (abs(z).lt.1e-8_dp) then
                        ispin = ispin + 1
                    endif
                else
                    ecore = z
                endif
                if ((ispin.eq.1).or.(ispin.eq.4)) then
                    ! aaaa/aa00
                    i = 2*i-1
                    if (ispin.eq.1) j = 2*j-1 ! so that it doesn't give -1
                    k = 2*k-1
                    if (ispin.eq.1) l = 2*l-1 ! so that it doesn't give -1
                elseif ((ispin.eq.2).or.(ispin.eq.5)) then
                    ! bbbb/bb00
                    i = 2*i
                    j = 2*j
                    k = 2*k
                    l = 2*l
                elseif (ispin.eq.3) then
                    ! aabb (chemical notation)
                    i = 2*i-1
                    k = 2*k-1
                    j = 2*j
                    l = 2*l
                endif
                if ((i.ne.0).and.(k.eq.0).and.(j.eq.0).and.(l.eq.0)) then
                    ! energies
                    arr(i) = z
                elseif ((i.ne.0).and.(k.ne.0).and.(j.eq.0).and.(l.eq.0)) then
                    ! <i|h|j>
                    tmat(i,k) = z
                    tmat(k,i) = z
                else
                    ! <ij|kl>
                    ! need to fill all 8 permutationally equivalent integrals
                    ! <ij|kl> = <ji|lk> = <kl|ij> = <lk|ji> = <kj|il>
                    ! = <li|jk> = <il|kj> = <jk|li>
                    umat(i,j,k,l) = z
                    umat(j,i,l,k) = z
                    umat(k,l,i,j) = z
                    umat(l,k,j,i) = z
                    umat(k,j,i,l) = z
                    umat(l,i,j,k) = z
                    umat(i,l,k,j) = z
                    umat(j,k,l,i) = z
                endif
            else
                if ((i.ne.0).and.(k.eq.0).and.(j.eq.0).and.(l.eq.0)) then
                    ! energies
                    arr(i) = z
                elseif ((i.ne.0).and.(k.ne.0).and.(j.eq.0).and.(l.eq.0)) then
                    ! <i|h|j>
                    tmat(i,k) = z
                    tmat(k,i) = z
                elseif ((i.eq.0).and.(k.eq.0).and.(j.eq.0).and.(l.eq.0)) then
                    ! core energy
                    ecore = z
                else
                    ! <ij|kl>
                    ! need to fill all 8 permutationally equivalent integrals
                    ! <ij|kl> = <ji|lk> = <kl|ij> = <lk|ji> = <kj|il>
                    ! = <li|jk> = <il|kj> = <jk|li>
                    umat(i,j,k,l) = z
                    umat(j,i,l,k) = z
                    umat(k,l,i,j) = z
                    umat(l,k,j,i) = z
                    umat(k,j,i,l) = z
                    umat(l,i,j,k) = z
                    umat(i,l,k,j) = z
                    umat(j,k,l,i) = z
                endif
            endif
        enddo
    ! complex integrals
    elseif (.not.trealinds) then
        ! if the integrals are in complex notation but still real
        ! the rotation is still done working only with real quantities
         complexint = .true.
         do
            read(8,*,end=199) cz,i,k,j,l
            if (abs(aimag(cz)).gt.1e-12_dp) then
                stop 'Error: Rotation cannot be done on complex orbitals'
            else
                complexint = .true.
            endif
            if (uhf.and.tmolpro) then
                ! Molpro writes out spatial orbitals indices which need 
                ! to be transfered to spin orbitals indices
                ! UHF FCIDUMP files in Molpro are written out as
                !1: aaaa
                !2: bbbb
                !3: aabb
                !4: aa
                !5: bb
                !with delimiters of 0.000000 0 0 0 0
                ! need to transform into spin orbital indices
                if ((i.eq.0).and.(k.eq.0).and.(j.eq.0).and.(l.eq.0).and.&
                    &(ispin.ne.6)) then
                    if (abs(cz).lt.1e-8_dp) then
                        ispin = ispin + 1
                    endif
                else
                    ecore = real(cz,dp)
                endif
                if ((ispin.eq.1).or.(ispin.eq.4)) then
                    ! aaaa/aa00
                    i = 2*i-1
                    if (ispin.eq.1) j = 2*j-1 ! so that it doesn't give -1
                    k = 2*k-1
                    if (ispin.eq.1) l = 2*l-1 ! so that it doesn't give -1
                elseif ((ispin.eq.2).or.(ispin.eq.5)) then
                    ! bbbb/bb00
                    i = 2*i
                    j = 2*j
                    k = 2*k
                    l = 2*l
                elseif (ispin.eq.3) then
                    ! aabb (chemical notation)
                    i = 2*i-1
                    k = 2*k-1
                    j = 2*j
                    l = 2*l
                endif
                if ((i.ne.0).and.(k.eq.0).and.(j.eq.0).and.(l.eq.0)) then
                    ! energies
                    arr(i) = real(cz,dp)
                elseif ((i.ne.0).and.(k.ne.0).and.(j.eq.0).and.(l.eq.0)) then
                    ! <i|h|j>
                    tmat(i,k) = real(cz,dp)
                    tmat(k,i) = real(cz,dp)
                else
                    ! <ij|kl>
                    ! need to fill all 4 permutationally equivalent integrals
                    ! <ij|kl> = <ji|lk> = (<kl|ij>)* = (<lk|ji>)* 
                    umat(i,j,k,l) = real(cz,dp)
                    umat(j,i,l,k) = real(cz,dp)
                    umat(k,l,i,j) = real(cz,dp)
                    umat(l,k,j,i) = real(cz,dp)
                endif
            else
                if ((i.ne.0).and.(k.eq.0).and.(j.eq.0).and.(l.eq.0)) then
                    ! energies
                    arr(i) = real(cz,dp)
                elseif ((i.ne.0).and.(k.ne.0).and.(j.eq.0).and.(l.eq.0)) then
                    ! <i|h|j>
                    tmat(i,k) = real(cz,dp)
                    tmat(k,i) = real(cz,dp)
                elseif ((i.eq.0).and.(k.eq.0).and.(j.eq.0).and.(l.eq.0)) then
                    ! core energy
                    ecore = real(cz,dp)
                else
                    ! <ij|kl>
                    ! need to fill all 4 permutationally equivalent integrals
                    ! <ij|kl> = <ji|lk> = (<kl|ij>)* = (<lk|ji>)* 
                    umat(i,j,k,l) = real(cz,dp)
                    umat(j,i,l,k) = real(cz,dp)
                    umat(k,l,i,j) = real(cz,dp)
                    umat(l,k,j,i) = real(cz,dp)
                endif
            endif
        enddo
    endif

    199 continue
    close(9)

    ! Molpro doesn't print out energies, so need to find them
    ! generate list which gives orbitals ordered according to 
    ! energy
    if (tmolpro) then
        call DiagFockElements(fdiag)
        arr(:) = fdiag(:)
        call OrderOrbs(energyorder)
    else
        call OrderOrbs(energyorder)
        call DiagFockElements(fdiag)
    endif

    ! Warning
    if (complexint) then
        write(6,*) '--------------------------------------------------------------'
        write(6,*) '******* Warning *******'
        write(6,*) 'Rotations is only performed using real parts of the integrals'
    endif


    ! Continue reading in sets of orbitals to rotate
    ! orbital indices refer to spatial orbitals if RHF and
    ! spin orbitals if UHF is considered
    open(9,file='rot_params',status='old',action='read')
    do l1=1,6
        read(9,*)
    enddo
    m = 0
    do
        ! each set of orbitals is specified by 2 rows
        ! if specifier=1: the second row is a list of orbital indices
        ! if specifier=2: the second row contains the start and end index 
        ! loc: specifies whether orbitals should be localised (if 1) or 
        ! delocalised (if 0) -> all values stored in localdelocal
        ! numb: number of orbitals in set to rotate -> all values
        ! stored in sets
        ! indices(:) : list of indices of orbitals to rotate
        ! startindex: start and of the set of orbitals
        ! to rotate -> indices of all orbitals to rotate stored in 
        ! rotate_list, rotate_list(m,:) contains all listed orbitals indices 
        ! (spatial orbitals for RHF and spin orbitals for UHF) for set m
        if (specifier.eq.1) then
            m = m + 1
            read(9,'(2I4)',end=299) loc,numb
            do l1=1,numb
                read(9,'(I4)',advance='no',end=299) indices(l1)
            enddo
            read(9,*,end=299)
            if (loc.eq.0) then
                localdelocal(m) = .true.
            else
                localdelocal(m) = .false.
            endif
            sets(m) = numb
            rotate_list(m,1:numb) = indices(1:numb)
        elseif (specifier.eq.2) then
            m = m + 1
            read(9,'(2I4)',end=299) loc,numb
            read(9,'(I4)',end=299) startindex
            if (loc.eq.1) then
                localdelocal(m) = .false.
            else
                localdelocal(m) = .true.
            endif
            sets(m) = numb
            do l1=1,numb
                rotate_list(m,l1) = startindex + ((l1-1)*1) 
            enddo
        endif
    enddo
    299 continue
    close(9)

    ! number of sets of orbitals
    m = m - 1
    nsets = m


    write(6,*) '---------------------------------------------------------------'
    write(6,*) 'Rotating the following sets of orbitals'
    ! index, (localise or delocalise?),number of orbitals in set, list of 
    ! obital indices
    write(6,*) 'Index of set, Localise?, Number of orbitals, Orbital Indices (spatial for RHF, spin for UHF)'
    do l1=1,nsets
        write(6,'(I5,L5,I7,A1)',advance='no') l1,localdelocal(l1),sets(l1),','
        do l2=1,sets(l1)
            write(6,'(I4)',advance='no') rotate_list(l1,l2)
        enddo
        write(6,*)
    enddo

    ! For UHF only orbitals with the same spin should be rotated amongst
    ! themselves
    if (uhf) then
        do l1=1,nsets
            do l2=2,sets(l1)
                if (mod(rotate_list(l1,(l2-1)),2).ne.mod(rotate_list(l1,l2),2)) then
                    stop 'Error in input options: only orbitals of the same spin&
                        & can be rotated amongst themselves'
                endif
            enddo
        enddo
    endif
    
    ! Generate a list of the rotations to perform
    ! Caution: Qchem orders the orbitals according to energy while
    ! Molpro orders them acordint to 
    ! the naximum number of pairs if norb*(norb-1) if 
    ! all orbitals were to be rotated in one set
    ! and 1 for the number of pairs in the set
    allocate(pairlist((norb*(norb-1)+1),nsets),stat=ierr)
    if (ierr.ne.0) then
        stop 'Error allocating list of pairs'
    endif
    pairlist(:,:) = 0

    ! if there are more then 2 orbitals in a set, the
    ! rotations have to be done self-consistenly, such
    ! that each iteraction corresponds to several pair-
    ! rotations such that each orbital is rotated with all
    ! of the reaining orbitals
    ! the only exception is when it is specified that degenerate
    ! orbitals should not be rotated amongst themselves, i.e.
    ! if trotdegen=.false.

    ! the requires lists of pairs are generated for each set
    write(6,*) '----------------------------------------------'
    write(6,*) 'Generating lists of pairs...'
    write(6,*) '----------------------------------------------'
    do l5=1,nsets
        call GeneratePairs(l5,sets(l5),rotate_list(l5,:),pairlist)
    enddo

    
    if (scheme.eq.1) then
        ! Rotate orbitals using selfinteraction criterion

        ! selfinteraction
        write(6,*) '----------------------------------------------------------'
        write(6,*) 'Selfinteractions before rotations:'
        write(6,*) 'Orbital, Selfinteraction'
        selfint(:) = 0.0_dp
        do l1=1,norb
            selfint(l1) = umat(l1,l1,l1,l1)
        enddo
        do l1=1,norb
            write(6,'(I3,3X,G25.12)') l1,selfint(l1)
        enddo
        selfint_curr = sum(selfint)
        write(6,*) 'Sum of selfinteractions:',selfint_curr
        write(6,*) '-----------------------------------------------------------'

        ! perform transformation
        transform(:,:) = 0.0_dp
        do l1=1,norb
            transform(l1,l1) = 1.0_dp
        enddo

 
        do l1=1,nsets
            write(6,*) 'Rotating Orbitals:',rotate_list(l1,1:sets(l1))
            write(6,'(A11,L3)') 'Localise ?',localdelocal(l1)
            rotatepairs = pairlist(1,l1)
            if (rotatepairs.eq.0) cycle
            iter = 0
            write(6,*) '--------------------------------------------------'
            write(6,*) 'Perfoming self-consistent rotation'
            do
                iter = iter + 1
                write(6,*) '----------------------------------------------------'
                write(6,*) 'Iteration:',iter
                selfint_prev = selfint_curr
                do l2=2,(2*rotatepairs),2
                    ! only 2 orbitals are rotated
                    transform(:,:) = 0.0_dp
                    do l3=1,norb
                        transform(l3,l3) = 1.0_dp
                    enddo
                    iumat = pairlist(l2,l1)
                    jumat = pairlist((l2+1),l1)
                    write(6,*) 'Rotating Pair:',iumat,jumat
                    call Rotate2Orbs(iumat,jumat,trans_2orbs_coeffs,selfint(iumat),&
                       selfint(jumat),localdelocal(l1))
                    transform(iumat,iumat) = trans_2orbs_coeffs(1,1)
                    transform(jumat,iumat) = trans_2orbs_coeffs(2,1)
                    transform(iumat,jumat) = trans_2orbs_coeffs(1,2)
                    transform(jumat,jumat) = trans_2orbs_coeffs(2,2)
                    
                    selfint_old = sum(selfint)
                    write(6,*) 'Sum of rotated self-interactions:',selfint_old
               
                    ! transform integrails
                    ! dgemm performs the operation
                    ! C = alpha * op(A) * op(B) + beta * op(c)
                    ! alpha,beta: scalars
                    ! op(a): m by k matrix
                    ! op(b): k by n matrix
                    ! op(C): m by n matrix
                    do i=1,norb
                        do j=1,norb
                            temp_inds(:,:) = 0.0_dp
                            call dgemm('T','N',norb,norb,norb,1.0_dp,transform(:,:),norb,&
                            &umat(1:norb,1:norb,i,j),norb,0.0_dp,&
                            &temp_inds(1:norb,1:norb),norb)
                            ! the transpose is needed
                            call dgemm('T','T',norb,norb,norb,1.0_dp,transform(:,:),norb,&
                            &temp_inds(1:norb,1:norb),norb,0.0_dp,&
                            &umat(1:norb,1:norb,i,j),norb)
                        enddo
                    enddo
                    do i=1,norb
                        do j=1,norb
                            temp_inds(:,:) = 0.0_dp
                            call dgemm('T','N',norb,norb,norb,1.0_dp,transform(:,:),norb,&
                            &umat(i,j,1:norb,1:norb),norb,0.0_dp,&
                            &temp_inds(1:norb,1:norb),norb)
                            ! the transpose is needed
                            call dgemm('T','T',norb,norb,norb,1.0_dp,transform(:,:),norb,&
                            &temp_inds(1:norb,1:norb),norb,0.0_dp,&
                            &umat(i,j,1:norb,1:norb),norb)
                        enddo
                    enddo

                    ! since the transformation matrix is overwritten each iteration
                    ! tmat needs to be transformed as well 
                    call dgemm('T','N',norb,norb,norb,1.0_dp,transform(:,:),norb,&
                    &tmat(1:norb,1:norb),norb,0.0_dp,&
                    &temp_tmat(1:norb,1:norb),norb)
                     call dgemm('T','T',norb,norb,norb,1.0_dp,transform(:,:),norb,&
                    &temp_tmat(1:norb,1:norb),norb,0.0_dp,&
                    &tmat(1:norb,1:norb),norb)

                    ! the same is true for fdiag
                    ! these 'transformed' orbital energy won't be correct 
                    ! because the fock matrix is not diagonal (unless it is the HF
                    ! basis) but they are sufficient for a guidance
                    call TransFDiag(fdiag)
                    arr(:) = fdiag(:)
                enddo

                selfint(:) = 0.0_dp
                do l2=1,norb
                    selfint(l2) = umat(l2,l2,l2,l2)
                enddo

                selfint_curr = sum(selfint)

                write(6,*) 'Selfinteractions of pervious and current iterations:'
                write(6,'(2G20.12)') selfint_prev,selfint_curr

                if (abs(selfint_curr-selfint_prev).lt.diff) then
                    if (localdelocal(l1)) then
                        if (selfint_curr.ge.selfint_prev) then
                            write(6,*) '-------------------------------------------'
                            write(6,*) 'Rotation is self-consistent !'
                            write(6,*) '-------------------------------------------'
                            exit 
                        else
                            cycle
                        endif
                    elseif (.not.localdelocal(l1)) then
                        if (selfint_curr.le.selfint_prev) then
                            write(6,*) '-------------------------------------------'
                            write(6,*) 'Rotation is self-consistent !'
                            write(6,*) '-------------------------------------------'
                            exit 
                        else
                            cycle
                        endif
                    endif
                endif
            enddo
        enddo

        ! selfinteraction
        write(6,*) '----------------------------------------------------------'
        write(6,*) 'Selfinteractions after rotations:'
        write(6,*) 'Orbital, Selfinteraction'
        selfint(:) = 0.0_dp
        do l1=1,norb
            selfint(l1) = umat(l1,l1,l1,l1)
        enddo
        do l1=1,norb
            write(6,'(I3,3X,G25.12)') l1,selfint(l1)
        enddo
        selfint_curr = sum(selfint)
        write(6,*) 'Sum of selfinteractions:',selfint_curr
        write(6,*) '-----------------------------------------------------------'

    elseif (scheme.eq.2) then
        ! Rotate orbitals using interactions \sum_{b in occ} <ib|ib|-<ib|bi>
        
        ! iteractions
        write(6,*) '----------------------------------------------------------'
        write(6,*) 'Interactions before rotations:'
        write(6,*) 'Orbital, Interaction'
        selfint(:) = 0.0_dp
        if (uhf) then
            do l1=1,norb
                do l2=1,nelec
                    l3 = energyorder(l2)
                    selfint(l1) = selfint(l1) +  umat(l1,l3,l1,l3)
                    if (mod(l1,2).eq.mod(l3,2)) then
                        selfint(l1) = selfint(l1) - umat(l1,l3,l3,l1)
                    endif
                enddo
            enddo
        else
            do l1=1,norb
                do l2=1,(nelec/2)
                    l3 = energyorder(l2)
                    selfint(l1) = selfint(l1) +  (2.0_dp*umat(l1,l3,l1,l3))&
                        & - umat(l1,l3,l3,l1)
                enddo
            enddo
        endif
        do l1=1,norb
            write(6,'(I3,3X,G25.12)') l1,selfint(l1)
        enddo
        selfint_curr = sum(selfint)
        write(6,*) 'Sum of Interactions:',selfint_curr
        write(6,*) '-----------------------------------------------------------'

        ! perform transformation
        transform(:,:) = 0.0_dp
        do l1=1,norb
            transform(l1,l1) = 1.0_dp
        enddo

        do l1=1,nsets
            write(6,*) 'Rotating Orbitals:',rotate_list(l1,1:sets(l1))
            write(6,'(A11,L3)') 'Maximise ?',localdelocal(l1)
            rotatepairs = pairlist(1,l1)
            if (rotatepairs.eq.0) cycle
            iter = 0
            write(6,*) '--------------------------------------------------'
            write(6,*) 'Perfoming self-consistent rotation'
            do
                iter = iter + 1
                write(6,*) '----------------------------------------------------'
                write(6,*) 'Iteration:',iter
                selfint_prev = selfint_curr
                do l2=2,(2*rotatepairs),2
                    ! only 2 orbitals are rotated
                    transform(:,:) = 0.0_dp
                    do l3=1,norb
                        transform(l3,l3) = 1.0_dp
                    enddo
                    iumat = pairlist(l2,l1)
                    jumat = pairlist((l2+1),l1)
                    write(6,*) 'Rotating Pair:',iumat,jumat
                    call RotateCoulombExchange(iumat,jumat,trans_2orbs_coeffs,selfint(iumat),&
                       selfint(jumat),localdelocal(l1))
                    transform(iumat,iumat) = trans_2orbs_coeffs(1,1)
                    transform(jumat,iumat) = trans_2orbs_coeffs(2,1)
                    transform(iumat,jumat) = trans_2orbs_coeffs(1,2)
                    transform(jumat,jumat) = trans_2orbs_coeffs(2,2)
                    
                    selfint_old = sum(selfint)
                    write(6,*) 'Sum of rotated Interactions:',selfint_old
               
                    ! transform integrails
                    ! dgemm performs the operation
                    ! C = alpha * op(A) * op(B) + beta * op(c)
                    ! alpha,beta: scalars
                    ! op(a): m by k matrix
                    ! op(b): k by n matrix
                    ! op(C): m by n matrix
                    do i=1,norb
                        do j=1,norb
                            temp_inds(:,:) = 0.0_dp
                            call dgemm('T','N',norb,norb,norb,1.0_dp,transform(:,:),norb,&
                            &umat(1:norb,1:norb,i,j),norb,0.0_dp,&
                            &temp_inds(1:norb,1:norb),norb)
                            ! the transpose is needed
                            call dgemm('T','T',norb,norb,norb,1.0_dp,transform(:,:),norb,&
                            &temp_inds(1:norb,1:norb),norb,0.0_dp,&
                            &umat(1:norb,1:norb,i,j),norb)
                        enddo
                    enddo
                    do i=1,norb
                        do j=1,norb
                            temp_inds(:,:) = 0.0_dp
                            call dgemm('T','N',norb,norb,norb,1.0_dp,transform(:,:),norb,&
                            &umat(i,j,1:norb,1:norb),norb,0.0_dp,&
                            &temp_inds(1:norb,1:norb),norb)
                            ! the transpose is needed
                            call dgemm('T','T',norb,norb,norb,1.0_dp,transform(:,:),norb,&
                            &temp_inds(1:norb,1:norb),norb,0.0_dp,&
                            &umat(i,j,1:norb,1:norb),norb)
                        enddo
                    enddo

                    ! since the transformation matrix is overwritten each iteration
                    ! tmat needs to be transformed as well 
                    call dgemm('T','N',norb,norb,norb,1.0_dp,transform(:,:),norb,&
                    &tmat(1:norb,1:norb),norb,0.0_dp,&
                    &temp_tmat(1:norb,1:norb),norb)
                     call dgemm('T','T',norb,norb,norb,1.0_dp,transform(:,:),norb,&
                    &temp_tmat(1:norb,1:norb),norb,0.0_dp,&
                    &tmat(1:norb,1:norb),norb)

                    ! the same is true for fdiag
                    ! these 'transformed' orbital energy won't be correct 
                    ! because the fock matrix is not diagonal (unless it is the HF
                    ! basis) but they are sufficient for a guidance
                    call TransFDiag(fdiag)
                    arr(:) = fdiag(:)
                enddo

                selfint(:) = 0.0_dp
                if (uhf) then
                    do l3=1,norb
                        do l2=1,nelec
                            l4 = energyorder(l2)
                            selfint(l3) = selfint(l3) +  umat(l3,l4,l3,l4)
                            if (mod(l3,2).eq.mod(l4,2)) then
                                selfint(l3) = selfint(l3) - umat(l3,l4,l4,l3)
                            endif
                        enddo
                    enddo
                else
                    do l3=1,norb
                        do l2=1,(nelec/2)
                            l4 = energyorder(l2)
                            selfint(l3) = selfint(l3) +  (2.0_dp*umat(l3,l4,l3,l4))&
                                & - umat(l3,l4,l4,l3)
                        enddo
                    enddo
                endif

                selfint_curr = sum(selfint)

                write(6,*) 'Interactions of pervious and current iterations:'
                write(6,'(2G20.12)') selfint_prev,selfint_curr

                if (abs(selfint_curr-selfint_prev).lt.diff) then
                    if (localdelocal(l1)) then
                        if (selfint_curr.ge.selfint_prev) then
                            write(6,*) '-------------------------------------------'
                            write(6,*) 'Rotation is self-consistent !'
                            write(6,*) '-------------------------------------------'
                            exit 
                        else
                            cycle
                        endif
                    elseif (.not.localdelocal(l1)) then
                        if (selfint_curr.le.selfint_prev) then
                            write(6,*) '-------------------------------------------'
                            write(6,*) 'Rotation is self-consistent !'
                            write(6,*) '-------------------------------------------'
                            exit 
                        else
                            cycle
                        endif
                    endif
                endif
            enddo
        enddo

        ! selfinteraction
        write(6,*) '----------------------------------------------------------'
        write(6,*) 'Interactions after rotations:'
        write(6,*) 'Orbital, Interaction'
        selfint(:) = 0.0_dp
        if (uhf) then
            do l3=1,norb
                do l2=1,nelec
                    l4 = energyorder(l2)
                    selfint(l3) = selfint(l3) +  umat(l3,l4,l3,l4)
                    if (mod(l3,2).eq.mod(l4,2)) then
                        selfint(l3) = selfint(l3) - umat(l3,l4,l4,l3)
                    endif
                enddo
            enddo
        else
            do l3=1,norb
                do l2=1,(nelec/2)
                    l4 = energyorder(l2)
                    selfint(l3) = selfint(l3) +  (2.0_dp*umat(l3,l4,l3,l4))&
                        & - umat(l3,l4,l4,l3)
                enddo
            enddo
        endif

        do l1=1,norb
            write(6,'(I3,3X,G25.12)') l1,selfint(l1)
        enddo
        selfint_curr = sum(selfint)
        write(6,*) 'Sum of Interactions:',selfint_curr
        write(6,*) '-----------------------------------------------------------'

    endif

    deallocate(temp_inds)

    deallocate(temp_tmat)

    ! writing out a new FCIUDMP_file file with rotated orbital integrals

    write(6,*) 'Writing FCIDUMP_rotated file...'

    open(12,FILE='FCIDUMP_rotated',status='unknown')
        
    write(12,'(2A6,I3,A7,I3,A5,I2,A)') '&FCI ','NORB=',norb,&
       ',NELEC=',nelec,',MS2=',ms2,','
    write(12,'(A9)',advance='no') 'ORBSYM='
    do i=1,norb
        write(12,'(I1,A1)',advance='no') 1,','
    enddo
    write(12,*) ''
    ! all orbital, k-point etx. symmetries are discarded
    ! since the orbitals have been rotated and the labels
    ! are thus no longer valid
    if (uhf) then
        write(12,'(A7,I1,A11)') 'ISYM=',1,' UHF=.TRUE.'
    else
        write(12,'(A7,I1,A12)') 'ISYM=',1,' UHF=.FALSE.'
    endif
    write(12,'(A5)') '&END'
    

    if (trealinds) then
        ! 2-electron integrals
        ! umat is in physcal notation, FCIDUMP is in
        ! chemical notation: <ij|kl> = (ik|jl)
        ! using permutational symmetry of real orbitals
        do i=1,norb
            do k=1,i
                do j=1,(i-1)
                    do l=1,j
                        if (abs(umat(i,j,k,l)).gt.1e-12_dp) then
                            write(12,'(1X,G20.12,4(I3))') umat(i,j,k,l),i,k,j,l
                        endif
                    enddo
                enddo
                j=i
                do l=1,k
                    if (abs(umat(i,j,k,l)).gt.1e-12_dp) then
                        write(12,'(1X,G20.12,4(I3))') umat(i,j,k,l),i,k,j,l
                    endif
                enddo
            enddo
        enddo

        ! 1-electron integrals
        do i=1,norb
            do j=1,norb
                if (abs(tmat(i,j)).gt.1e-12_dp) then
                    write(12,'(1X,G20.12,4(I3))') tmat(i,j),i,j,0,0
                endif
            enddo
        enddo

        ! energies
        do i=1,norb
            write(12,'(1X,G20.12,4(I3))') arr(i),i,0,0,0
        enddo

        ! core energy
        write(12,'(1X,G20.12,4(I3))') ecore,0,0,0,0

    elseif (.not.trealinds) then
        ! using permutational symmetry of complex orbitals
        do i=1,norb
            do j=1,i
                do k=1,i
                    do l=1,i
                        if ((k.gt.i).and.(l.gt.j)) cycle
                        if (abs(umat(i,j,k,l)).gt.1e-12_dp) then
                            write(12,*) cmplx(umat(i,j,k,l),0.0_dp,dp),i,k,j,l
                        endif
                    enddo
                enddo
            enddo
        enddo
 
        ! 1-electron integrals
        do i=1,norb
            do j=1,norb
                if (abs(tmat(i,j)).gt.1e-12_dp) then
                    write(12,*) cmplx(tmat(i,j),0.0_dp,dp),i,j,0,0
                endif
            enddo
        enddo

        ! energies
        do i=1,norb
            write(12,*) cmplx(arr(i),0.0_dp,dp),i,0,0,0
        enddo

        ! core energy
        write(12,*) cmplx(ecore,0.0_dp,dp),0,0,0,0

    endif

    close(12)

    deallocate(localdelocal)
    deallocate(sets)
    deallocate(rotate_list)
    deallocate(pairlist)
    deallocate(umat)
    deallocate(tmat)
    deallocate(arr)
    deallocate(selfint)
    deallocate(transform)
    deallocate(trans_2orbs_coeffs)
 
contains

    subroutine Rotate2Orbs(i,j,trans_2orbs_coeffs,selfintorb1,selfintorb2,localdelocal)

        ! This routine takes two orbitals i,j, and rotates them in order to maximally localise these
        ! It employs an Edminston-Ruedenberg type localisation which maximises the self-interaction
        ! \sum_{i=1}^{2} \sum_{r,s,u,v} (c_{ir})*(c_{is})*c_{iu}c_{iv} <p_{i}p_{i}|u|p_{i}p_{i}>
        ! where p_{i} are the original NOs
        ! The the coefficients c are given by the following matrix:
        ! c =  cos a   sin a
        !      -sin a  cos a
        ! Then angle a is found by differentiating and setting it equal to 0 which gives
        ! the following analytical expression of the form
        ! tan a = -x/y
        ! where x and y are sums of the original NO four index inegrals
        real(dp), allocatable, intent(inout) :: trans_2orbs_coeffs(:,:)
        real(dp), intent(inout) :: selfintorb1,selfintorb2
        real(dp) :: alpha2(17)
        !real(dp) :: secondderiv(2)
        real(dp) :: selfinteractions(17)
        real(dp) :: coeffcos,coeffsin,maxint
        integer :: maxangle(1)
        integer :: indicesij(2)
        integer, intent(in) :: i,j
        integer :: l1,l2,l3,l4,l5
        logical, intent(in) :: localdelocal

        indicesij(1) = i
        indicesij(2) = j
        trans_2orbs_coeffs(:,:) = 0.0_dp

        ! Umat(i,j,k,l) contains the four-index integrals
        ! <ij|kl> (physical notation) in the NO basis

        coeffcos = Umat(i,i,i,j) + Umat(i,i,j,i) + Umat(i,j,i,i) &
            & - Umat(i,j,j,j) + Umat(j,i,i,i) - Umat(j,i,j,j) &
            & - Umat(j,j,i,j) - Umat(j,j,j,i)

        coeffsin = -Umat(i,i,i,i) + Umat(i,i,j,j) + Umat(i,j,i,j) &
            & + Umat(i,j,j,i) + Umat(j,i,i,j) + Umat(j,i,j,i) &
            & + Umat(j,j,i,i) - Umat(j,j,j,j)

        ! atan return a value in [-pi/2,pi/2]
        ! because of the 4*alpha in the equation there are 8 distinct solutions
        ! i.e. in the range 0,2*pi
        ! i.e. possible solutions are separated by (2*pi/8)=pi/4
        ! for safety 16 solutions are evaluated
        alpha2(9) = atan((-coeffcos/coeffsin))
        alpha2(9) = alpha2(9)/4.0_dp
        do l1=8,1,-1
            alpha2(l1) = alpha2(l1+1) - (pi/4.0_dp)
        enddo
        do l1=10,17
            alpha2(l1) = alpha2(l1-1) + (pi/4.0_dp)
        enddo
        
        !! second derivatives to find maximum (necessary since the minimum, i.e. fully delocalised
        !! orbitals satisfy the same conditions
        !secondderiv(1) = (4.0_dp*coeffsin*cos(4.0_dp*alpha(1))) - (4.0_dp*coeffcos*sin(4.0_dp*alpha(1)))
        !secondderiv(2) = (4.0_dp*coeffsin*cos(4.0_dp*alpha(2))) - (4.0_dp*coeffcos*sin(4.0_dp*alpha(2)))

        ! compute selfinteractions to check which one is largest
        ! this is a better measure than the second derivatives
        selfinteractions(:) = 0.0_dp 

        do l1=1,17
            trans_2orbs_coeffs(1,1) = cos(alpha2(l1))
            trans_2orbs_coeffs(2,1) = sin(alpha2(l1))
            trans_2orbs_coeffs(1,2) = -sin(alpha2(l1))
            trans_2orbs_coeffs(2,2) = cos(alpha2(l1))

            do l2=1,2
                do l3=1,2
                    do l4=1,2
                        do l5=1,2
                            selfinteractions(l1) = selfinteractions(l1) + trans_2orbs_coeffs(l2,1)&
                                &*trans_2orbs_coeffs(l3,1)*&
                                &trans_2orbs_coeffs(l4,1)*trans_2orbs_coeffs(l5,1)*&
                                &Umat(indicesij(l2),indicesij(l3),indicesij(l4),indicesij(l5))
                        enddo
                    enddo
                enddo
            enddo
            do l2=1,2
                do l3=1,2
                    do l4=1,2
                        do l5=1,2
                            selfinteractions(l1) = selfinteractions(l1) + trans_2orbs_coeffs(l2,2)&
                                &*trans_2orbs_coeffs(l3,2)*&
                                &trans_2orbs_coeffs(l4,2)*trans_2orbs_coeffs(l5,2)*&
                                &Umat(indicesij(l2),indicesij(l3),indicesij(l4),indicesij(l5))
                        enddo
                    enddo
                enddo
            enddo
        enddo

        !do l1=1,17
        !    write(6,'(I3,1X,5(G20.12))') l1,alpha2(l1),tan(alpha2(l1)*4.0_dp),cos(alpha2(l1)),&
        !        &sin(alpha2(l1)),selfinteractions(l1)
        !enddo

        ! choose the angle which maximises the selfinteractions
        if (.not.localdelocal) then
            ! maximally delocalised
            maxangle = minloc(selfinteractions)
            maxint = minval(selfinteractions)
        elseif (localdelocal) then
            ! maximally localised
            maxangle = maxloc(selfinteractions)
            maxint = maxval(selfinteractions)
        endif


        ! return transformatin coefficients
        trans_2orbs_coeffs(1,1) = cos(alpha2(maxangle(1)))
        trans_2orbs_coeffs(2,1) = sin(alpha2(maxangle(1)))
        trans_2orbs_coeffs(1,2) = -sin(alpha2(maxangle(1)))
        trans_2orbs_coeffs(2,2) = cos(alpha2(maxangle(1)))
 
        ! new sefl-interactions for transformed orbitals
        selfintorb1 = 0.0_dp
        selfintorb2 = 0.0_dp
        do l2=1,2
            do l3=1,2
                do l4=1,2
                    do l5=1,2
                        selfintorb1 = selfintorb1 + trans_2orbs_coeffs(l2,1)&
                            &*trans_2orbs_coeffs(l3,1)*&
                            &trans_2orbs_coeffs(l4,1)*trans_2orbs_coeffs(l5,1)*&
                            &Umat(indicesij(l2),indicesij(l3),indicesij(l4),indicesij(l5))
                        selfintorb2 = selfintorb2 + trans_2orbs_coeffs(l2,2)&
                            &*trans_2orbs_coeffs(l3,2)*&
                            &trans_2orbs_coeffs(l4,2)*trans_2orbs_coeffs(l5,2)*&
                            &Umat(indicesij(l2),indicesij(l3),indicesij(l4),indicesij(l5))
                    enddo
                enddo
            enddo
        enddo


    endsubroutine Rotate2Orbs


    subroutine GeneratePairs(numberset,numb,orbitalsset,pairs)

        integer, intent(in) :: numb,numberset
        integer, intent(in) :: orbitalsset(:)
        integer, allocatable, intent(inout) :: pairs(:,:)
    !    integer :: frac!,shift
        integer :: l1,l2,npairs,m!,n,p

        ! This routine generates the list of pairs which is needed for performing
        ! the rotations self-consistenly by step-wise 2-orbital rotations

    !    if (mod(frac,2).eq.0) then
    !        frac = numb - 1
    !    elseif (mod(frac,2).eq.1) then
    !        frac = numb - 1
    !    endif

        npairs = 0

        if (numb.eq.2) then
            ! if there are 2 orbitals in the set only one rotation is needed
            npairs = npairs + 1
            pairs((1+1),numberset) = orbitalsset(1)
            pairs((2+1),numberset) = orbitalsset(2)
        else
            ! if there are mor orbitals in the set the rotations have to
            ! be done self-consistently where each iteraction consists of 
            ! a set of rotations such that each orbital has been rotated
            ! with all of the remainder
            ! the only exception to this is when it is specified that set
            ! degenerate orbitals should not be rotated amongst themselfes
            ! i.e. if trotdegen-.false.
            m = 0
            do l1=1,numb
                do l2=(l1+1),numb
                    m = m + 1
                    npairs = npairs + 1
                    pairs((m*2),numberset) = orbitalsset(l1)
                    pairs(((m*2)+1),numberset) = orbitalsset(l2)
                    if ((.not.trotdegen).and.(abs(arr(orbitalsset(l1))-arr(orbitalsset(l2)))&
                        &.lt.1e-12_dp)) then
                        m = m - 1
                        npairs = npairs - 1
                    endif
                enddo
            enddo
      !      m = 0
      !      p = 0
      !      do l1=1,numb!1,numb
      !          !shift = frac*(l1-1)
      !          m = m + 1
      !          n = 0
      !          p = 0
      !          do l2=(l1+1),numb
      !              !m = m + 1
      !              n = n + 1
      !              npairs = npairs + 1
      !              shift = frac*(n-1)
      !              if (n.gt.2) then 
      !                  p = p + 1
      !                  shift = shift - p
      !              endif
      !              write(6,*) l1,orbitalsset(l1),l2,orbitalsset(l2),shift,((shift+((m)))*2)
      !              pairs(((shift+m)*2),numberset) = orbitalsset(l1)
      !              pairs((((shift+m)*2)+1),numberset) = orbitalsset(l2)
!     !                pairs(((2*m)+(shift)),numberset) = orbitalsset(l1)
!     !               pairs(((2*m)+(shift)+1),numberset) = orbitalsset(l2)
      !          enddo
      !      enddo
        endif

        pairs(1,numberset) = npairs



    endsubroutine GeneratePairs

    subroutine OrderOrbs(energyorder)


        ! This subroutine orders the orbitals 
        ! It returns a list of the orbitals ordered according to energy
        
        integer, allocatable, intent(inout) :: energyorder(:)
        real(dp), allocatable :: sortedarr(:)
        integer :: l1,tempind
        real(dp) :: temp
        logical :: swapped

        allocate(sortedarr(norb))

        sortedarr(:) = arr(:)
        
        do l1=1,norb
            energyorder(l1) = l1
        enddo
        
        ! sorting according to ascending energy
        do
            swapped=.false.
            do l1=2,norb
                if (arr(l1).lt.arr(l1-1)) then
                    temp = arr(l1-1)
                    tempind = energyorder(l1-1)
                    arr(l1-1) = arr(l1)
                    energyorder(l1-1) = energyorder(l1)
                    arr(l1) = temp
                    energyorder(l1) = tempind
                    swapped=.true.
                endif
            enddo
            if (.not.swapped) exit
        enddo

        deallocate(sortedarr)


    endsubroutine OrderOrbs

    subroutine DiagFockElements(fdiag)

        real(dp), allocatable,intent(inout) :: fdiag(:)
        integer :: i,j,b,g

        ! This routine calculates the fock matrix from the initial read in values
        ! of umat and tmat
        ! in the HF orbital basis the fock matrix is diagonal with 
        ! with diagonal elements being 
        ! e_i = f__ii = h_ii + \sum_[b} <ib|ib> - <ib|bi>

        ! Only the diagonal Fock matrix elments are calculated since only 
        ! these are used: in an HF basis the Fock matrix is diagonal and so
        ! this will be exact, for other bases this will not be correct but
        ! it is sufficient to give a guideline for the orbital energies
        ! of the rotated orbitals

        if (tmolpro) then
            if (uhf) then
                ! nelec orbitals with lowest indices will always be 
                ! the occupied orbitals
                fdiag(:) = 0.0_dp
                ! the one-electron integrals 
                do i=1,norb
                    fdiag(i) = tmat(i,i)
                enddo

                ! the two-electron integrals
                ! coulomb contribution
                do i=1,norb
                    !j = energyorder(i)
                    do b=1,nelec
                        !g = energyorder(b)
                        fdiag(i) = fdiag(i) + umat(i,b,i,b)
                    enddo
                enddo

                ! exchange contribution 
                ! this is zero of electrons are of different spin
                ! all odd integer spin orbitals are of alpha spin
                ! all even integer spin orbitals are of beta spin
                do i=1,norb
                    !j = energyorder(i)
                    do b=1,nelec
                        !g = energyorder(j)
                        if ((mod(i,2).eq.mod(b,2))) then
                            fdiag(i) = fdiag(i) - umat(i,b,b,i)
                        endif
                    enddo
                enddo
            else
                ! nelec orbitals with lowest indices will always be 
                ! the occupied orbitals
                fdiag(:) = 0.0_dp
                ! the one-electron integrals 
                do i=1,norb
                    fdiag(i) = tmat(i,i)
                enddo

                ! the two-electron integrals
                ! coulomb contribution
                do i=1,norb
                    !j = energyorder(i)
                    do b=1,(nelec/2)
                        !g = energyorder(b)
                        fdiag(i) = fdiag(i) + (2.0_dp*umat(i,b,i,b))
                    enddo
                enddo

                ! exchange contribution 
                ! this is zero of electrons are of different spin
                ! all odd integer spin orbitals are of alpha spin
                ! all even integer spin orbitals are of beta spin
                do i=1,norb
                    !j = energyorder(i)
                    do b=1,(nelec/2)
                        !g = energyorder(j)
                        fdiag(i) = fdiag(i) - umat(i,b,b,i)
                    enddo
                enddo
 
            endif
        else
            if (uhf) then
                fdiag(:) = 0.0_dp
                ! the one-electron integrals 
                do i=1,norb
                    fdiag(i) = tmat(i,i)
                enddo

                ! the two-electron integrals
                ! coulomb contribution
                do i=1,norb
                    j = energyorder(i)
                    do b=1,nelec
                        g = energyorder(b)
                        fdiag(j) = fdiag(j) + umat(j,g,j,g)
                    enddo
                enddo

                ! exchange contribution 
                ! this is zero of electrons are of different spin
                ! all odd integer spin orbitals are of alpha spin
                ! all even integer spin orbitals are of beta spin
                do i=1,norb
                    j = energyorder(i)
                    do b=1,nelec
                        g = energyorder(b)
                        if ((mod(j,2).eq.mod(g,2))) then
                            fdiag(j) = fdiag(j) - umat(j,g,g,j)
                        endif
                    enddo
                enddo
            else
                fdiag(:) = 0.0_dp
                ! the one-electron integrals 
                do i=1,norb
                    fdiag(i) = tmat(i,i)
                enddo

                ! the two-electron integrals
                ! coulomb contribution
                do i=1,norb
                    j = energyorder(i)
                    do b=1,(nelec/2)
                        g = energyorder(b)
                        fdiag(j) = fdiag(j) + (2.0_dp*umat(j,g,j,g))
                    enddo
                enddo

                ! exchange contribution 
                ! this is zero of electrons are of different spin
                ! all odd integer spin orbitals are of alpha spin
                ! all even integer spin orbitals are of beta spin
                do i=1,norb
                    j = energyorder(i)
                    do b=1,(nelec/2)
                        g = energyorder(b)
                        fdiag(j) = fdiag(j) - umat(j,g,g,j)
                    enddo
                enddo
            endif
        endif

    endsubroutine DiagFockElements

    subroutine TransFDiag(fdiag)

        ! This subroutine transforms the diagonal elements of the fock
        ! matrix to give an estimate of the orbital energies 

        real(dp), allocatable, intent(inout) :: fdiag(:)
        real(dp), allocatable :: tempfdiag(:)
        integer :: i,j

        allocate(tempfdiag(norb))

        tempfdiag(:) = fdiag(:)
        fdiag(:) = 0.0_dp

        do i=1,norb
            do j=1,norb
                fdiag(i) = fdiag(i) + transform(j,i)*tempfdiag(j)&
                    &*transform(j,i)
            enddo
        enddo

        deallocate(tempfdiag)


    endsubroutine TransFDiag

    subroutine RotateCoulombExchange(i,j,trans_2orbs_coeffs,selfintorb1,selfintorb2,localdelocal)

        ! This routine takes two orbitals i,j, and rotates them in order to maximally localise these
        ! It employs a localisation scheme which maximises/minimises the electron interaction
        ! terms  in the HF equation,i.e. the combined coulomb and exchange interactions
        ! \sum_{i=1}^{2} \sum_{r,s} \sum_{j in occ} (c_{ir})*c_{is} <p_{i}p_{j}|u|p_{i}p_{j}>) 
        ! - (c_{ir})*c_{is} <p_{i}p_{j}|u|p_{i}p_{j}>)
        ! where p_{i} are the original NOs
        ! The the coefficients c are given by the following matrix:
        ! c =  cos a   sin a
        !      -sin a  cos a
        ! Then angle a is found by evaluating the quantity for a discrete set of angles
        ! sampled at a high enough frequency (Nyqusist's theorem: twice the highest
        ! fourier component frequency)

        real(dp), allocatable, intent(inout) :: trans_2orbs_coeffs(:,:)
        real(dp), intent(inout) :: selfintorb1,selfintorb2
        real(dp) :: alpha(32)
        real(dp) :: selfinteractions(32)
        real(dp) :: maxint
        integer :: maxangle(1)
        integer :: indicesij(2)
        integer, intent(in) :: i,j
        integer :: l1,l2,l3,l4,l5
        logical, intent(in) :: localdelocal

        indicesij(1) = i
        indicesij(2) = j
        trans_2orbs_coeffs(:,:) = 0.0_dp


        ! the angle alpha is sampled in the region 0,4*pi
        alpha(1) = 0
        do l1=2,32
            alpha(l1) = alpha(l1-1) + (2.0_dp*pi/32.0_dp)
        enddo


        ! compute coulomb and exchange interactions 
        ! to find the maximum is largest
        selfinteractions(:) = 0.0_dp 

        do l1=1,32
            trans_2orbs_coeffs(1,1) = cos(alpha(l1))
            trans_2orbs_coeffs(2,1) = sin(alpha(l1))
            trans_2orbs_coeffs(1,2) = -sin(alpha(l1))
            trans_2orbs_coeffs(2,2) = cos(alpha(l1))

            if (uhf) then

                ! first orbital
                do l2=1,2
                    do l3=1,2
                        do l4=1,nelec
                            ! for correct reference determinant
                            l5 = energyorder(l4)
                            ! coulomb integrals
                            selfinteractions(l1) = selfinteractions(l1) + (trans_2orbs_coeffs(l2,1)&
                                &*trans_2orbs_coeffs(l3,1)*&
                                &Umat(indicesij(l2),l5,indicesij(l3),l5))
                            ! exchange integrals
                            ! only non-zero for electrons with the same spin
                            if ((mod(indicesij(l3),2).eq.mod(indicesij(l2),2)).and.&
                                &(mod(indicesij(l3),2).eq.mod(l5,2))) then
                                selfinteractions(l1) = selfinteractions(l1) - (trans_2orbs_coeffs(l2,1)&
                                    &*trans_2orbs_coeffs(l3,1)*&
                                    &Umat(indicesij(l2),l5,l5,indicesij(l3)))
                            endif
                        enddo
                    enddo
                enddo
                ! the sum of both orbitals is invariant
                ! so only one orbital can be maximised/minimised 
                ! while the other one will be minimised/maximised
    !            ! second orbital
    !            do l2=1,2
    !                do l3=1,2
    !                    do l4=1,nelec
    !                        ! correct reference determinant
    !                        l5 = energyorder(l4)
    !                        ! coulomb integrals
    !                        selfinteractions(l1) = selfinteractions(l1) + (trans_2orbs_coeffs(l2,2)&
    !                            &*trans_2orbs_coeffs(l3,2)*&
    !                            &Umat(indicesij(l2),l5,indicesij(l3),l5))
    !                        ! eschange integrals
    !                        ! only non-zero for electrons with the same spin
    !                        if ((mod(indicesij(l3),2).eq.mod(indicesij(l2),2)).and.&
    !                            &(mod(indicesij(l3),2).eq.mod(l5,2))) then
    !                            selfinteractions(l1) = selfinteractions(l1) - (trans_2orbs_coeffs(l2,2)&
    !                                &*trans_2orbs_coeffs(l3,2)*&
    !                                &Umat(indicesij(l2),l5,l5,indicesij(l3)))
    !                        endif
    !                    enddo
    !                enddo
    !            enddo
            else
                ! first orbital
                do l2=1,2
                    do l3=1,2
                        do l4=1,(nelec/2)
                            ! for correct reference determinant
                            l5 = energyorder(l4)
                            ! coulomb integrals
                            selfinteractions(l1) = selfinteractions(l1) + (trans_2orbs_coeffs(l2,1)&
                                &*trans_2orbs_coeffs(l3,1)*2.0_dp*&
                                &Umat(indicesij(l2),l5,indicesij(l3),l5))
                            ! exchange integrals
                            ! only non-zero for electrons with the same spin
                            selfinteractions(l1) = selfinteractions(l1) - (trans_2orbs_coeffs(l2,1)&
                                &*trans_2orbs_coeffs(l3,1)*&
                                &Umat(indicesij(l2),l5,l5,indicesij(l3)))
                        enddo
                    enddo
                enddo
                ! the sum of both orbitals is invariant
                ! so only one orbital can be maximised/minimised 
                ! while the other one will be minimised/maximised
    !            ! second orbital
    !            do l2=1,2
    !                do l3=1,2
    !                    do l4=1,(nelec/2)
    !                        ! correct reference determinant
    !                        l5 = energyorder(l4)
    !                        ! coulomb integrals
    !                        selfinteractions(l1) = selfinteractions(l1) + (trans_2orbs_coeffs(l2,2)&
    !                            &*trans_2orbs_coeffs(l3,2)*2.0_dp*&
    !                            &Umat(indicesij(l2),l5,indicesij(l3),l5))
    !                        ! eschange integrals
    !                        ! only non-zero for electrons with the same spin
    !                        selfinteractions(l1) = selfinteractions(l1) - (trans_2orbs_coeffs(l2,2)&
    !                            &*trans_2orbs_coeffs(l3,2)*&
    !                            &Umat(indicesij(l2),l5,l5,indicesij(l3)))
    !                    enddo
    !                enddo
    !            enddo
            endif
        enddo

        !do l1=1,32
        !    write(6,'(I3,1X,5(G20.12))') l1,alpha(l1),cos(alpha(l1)),sin(alpha(l1)),&
        !        &selfinteractions(l1)
        !enddo

        ! choose the angle which maximises the selfinteractions
        if (.not.localdelocal) then
            ! maximum value
            maxangle = minloc(selfinteractions)
            maxint = minval(selfinteractions)
        elseif (localdelocal) then
            ! minimum values
            maxangle = maxloc(selfinteractions)
            maxint = maxval(selfinteractions)
        endif


        ! return transformatin coefficients
        trans_2orbs_coeffs(1,1) = cos(alpha(maxangle(1)))
        trans_2orbs_coeffs(2,1) = sin(alpha(maxangle(1)))
        trans_2orbs_coeffs(1,2) = -sin(alpha(maxangle(1)))
        trans_2orbs_coeffs(2,2) = cos(alpha(maxangle(1)))
 
 
        ! quantity \sum_{k in occ} <ik|ik> - <ik|ki>
        ! for rotated orbitals
        selfintorb1 = 0.0_dp
        selfintorb2 = 0.0_dp
        if (uhf) then

            ! first orbital
            do l2=1,2
                do l3=1,2
                    do l4=1,nelec
                        ! for correct reference determinant
                        l5 = energyorder(l4)
                        ! coulomb integrals
                        selfintorb1 = selfintorb1 + (trans_2orbs_coeffs(l2,1)&
                            &*trans_2orbs_coeffs(l3,1)*&
                            &Umat(indicesij(l2),l5,indicesij(l3),l5))
                        ! exchange integrals
                        ! only non-zero for electrons with the same spin
                        if ((mod(indicesij(l3),2).eq.mod(indicesij(l2),2)).and.&
                            &(mod(indicesij(l3),2).eq.mod(l5,2))) then
                            selfintorb1 = selfintorb1 - (trans_2orbs_coeffs(l2,1)&
                                &*trans_2orbs_coeffs(l3,1)*&
                                &Umat(indicesij(l2),l5,l5,indicesij(l3)))
                        endif
                    enddo
                enddo
            enddo
            ! second orbital
            do l2=1,2
                do l3=1,2
                    do l4=1,nelec
                        ! correct reference determinant
                        l5 = energyorder(l4)
                        ! coulomb integrals
                        selfintorb2 = selfintorb2 + (trans_2orbs_coeffs(l2,2)&
                            &*trans_2orbs_coeffs(l3,2)*&
                            &Umat(indicesij(l2),l5,indicesij(l3),l5))
                        ! eschange integrals
                        ! only non-zero for electrons with the same spin
                        if ((mod(indicesij(l3),2).eq.mod(indicesij(l2),2)).and.&
                            &(mod(indicesij(l3),2).eq.mod(l5,2))) then
                            selfintorb2 = selfintorb2 - (trans_2orbs_coeffs(l2,2)&
                                &*trans_2orbs_coeffs(l3,2)*&
                                &Umat(indicesij(l2),l5,l5,indicesij(l3)))
                        endif
                    enddo
                enddo
            enddo
        else
            ! first orbital
            do l2=1,2
                do l3=1,2
                    do l4=1,(nelec/2)
                        ! for correct reference determinant
                        l5 = energyorder(l4)
                        ! coulomb integrals
                        selfintorb1 = selfintorb1 + (trans_2orbs_coeffs(l2,1)&
                            &*trans_2orbs_coeffs(l3,1)*2.0_dp*&
                            &Umat(indicesij(l2),l5,indicesij(l3),l5))
                        ! exchange integrals
                        ! only non-zero for electrons with the same spin
                        selfintorb1 = selfintorb1 - (trans_2orbs_coeffs(l2,1)&
                            &*trans_2orbs_coeffs(l3,1)*&
                            &Umat(indicesij(l2),l5,l5,indicesij(l3)))
                    enddo
                enddo
            enddo
            ! second orbital
            do l2=1,2
                do l3=1,2
                    do l4=1,(nelec/2)
                        ! correct reference determinant
                        l5 = energyorder(l4)
                        ! coulomb integrals
                        selfintorb2 = selfintorb2 + (trans_2orbs_coeffs(l2,2)&
                            &*trans_2orbs_coeffs(l3,2)*2.0_dp*&
                            &Umat(indicesij(l2),l5,indicesij(l3),l5))
                        ! eschange integrals
                        ! only non-zero for electrons with the same spin
                        selfintorb2 = selfintorb2 - (trans_2orbs_coeffs(l2,2)&
                            &*trans_2orbs_coeffs(l3,2)*&
                            &Umat(indicesij(l2),l5,l5,indicesij(l3)))
                    enddo
                enddo
            enddo
        endif


    endsubroutine RotateCoulombExchange

end program main
