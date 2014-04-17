program main
! This program reads in an FCIDMP input file and performs rotation of the 
! orbitals in order to break symmetry and maximally localise or delocalise 
! the rotated orbitals.
! The rotations are performed pair-wise by calling the Rotate2Orbs routine.
! compile with: f95 -F2008 -o fcidump_instability.x fcidump_instability.f90 -llapack -lblas -latlas
! / -lacml

implicit none

! constant date

!integer, parameter :: sp = selected_real_kind(6,37)
integer, parameter :: dp = selected_real_kind(15,307)
!integer, parameter :: qp = selected_real_kind(33,4931)
!integer, parameter :: int32 = selected_int_kind(8)
integer, parameter :: int64 = selected_int_kind(15)

real(dp), parameter :: pi = 3.1415926535897932384626433832795028841971693993751_dp
complex(dp), parameter :: imaginaryi = (0.0_dp,1.0_dp)

integer :: norb,nelec,ms2,isym,syml(1000),iuhf,symlz(1000),nprop(3)
integer :: propbitlen
logical :: uhf
real(dp), allocatable :: umat(:,:,:,:),tmat(:,:),arr(:),temp_inds(:,:)
real(dp), allocatable :: temp_tmat(:,:)
real(dp), allocatable :: asinglet(:,:),atriplet(:,:),bsinglet(:,:),btriplet(:,:)
real(dp), allocatable :: xsinglet(:,:),xtriplet(:,:)
real(dp), allocatable :: xsinglet_breaktimerev(:,:),xtriplet_breaktimerev(:,:)
real(dp), allocatable :: xsinglet_eval(:),xtriplet_eval(:)
real(dp), allocatable :: xsinglet_breaktimerev_eval(:),xtriplet_breaktimerev_eval(:)
real(dp), allocatable :: smat(:,:),temp_smat(:,:),temp2_smat(:,:)
real(dp), allocatable :: scr(:),smat_eval(:)
real(dp) :: ecore
complex(dp) :: cz
complex(dp), allocatable :: coeffs(:,:),temptrans(:,:)
complex(dp), allocatable :: grade1(:,:),grade2(:,:,:,:)
real(dp) :: z
integer :: loc(2)
real(dp) :: ratio,angle,hf_energy
real(dp), allocatable :: inst_all(:,:),transarray(:)
integer :: nmonoex,scheme
integer(int64) :: orbsym(1000)
integer :: i,j,k,l,a,ierr,ispin,m,lscr
integer :: l1,l2,l3
real(dp), allocatable :: transform(:,:),fdiag(:)
integer, allocatable :: energyorder(:)
integer, allocatable :: monoexcitations(:,:)
logical :: exists,trealinds,tmolpro,complexint
namelist /fci/ norb,nelec,ms2,orbsym,isym,iuhf,uhf,syml,symlz,propbitlen,nprop
    
    ! Defaults for rot_params input file
    trealinds = .true.
    tmolpro = .false.

    ! Remaining defaults 
    uhf = .false.
    propbitlen = 0
    nprop = 0
    iuhf = 0

    write(6,*) '---------------------------------------------------'
    write(6,*) 'Calculating instability of solution in FCIDUMP file'
    write(6,*) '---------------------------------------------------'
 
    ! Read in input options for rotation
    inquire(file='inst_params',exist=exists)
    if (.not.exists) then
        stop 'inst_params input file does not exist'
    endif
    write(6,*) 'Reading in inst_params input file...'
    open(9,file='inst_params',status='old',action='read')
    read(9,*) trealinds      ! integrals in FCIDUMP file are real (if .true.) and complex 
                             ! (if .false.)
    read(9,*) tmolpro        ! Molpro FCIDUMP file (if .true.) or Qchem FCIDUMP file 
                             ! (if .false.)
                             ! each other (if.true.) and otherwise not (if.false.)
    read(9,*) scheme         ! scheme which is used for rotating the orbitals in order to find the lower energy solution, if scheme=1 a pairwise orbital rotation is performed, if scheme=2: the transformation is performed using the eigenvector followed by a symmetric Loewdin orthogonalisation, if scheme=3 a further RHF calculation is performed 
    close(9)


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

 
    ! in order to find any possible instablity in the RHF/UHF solutions
    ! all possible monoexcitations have to be considered

    ! number of all possible monexcitations
    if (uhf) then
        nmonoex = nelec*(norb-nelec)
    elseif (.not.uhf) then
        nmonoex = (nelec/2)*(norb-(nelec/2))
    endif
   
    allocate(monoexcitations(nmonoex,2),stat=ierr)
    if (ierr.ne.0) then
        stop 'Error allocating list for monoexcitations'
    endif

    write(6,*) '------------------------------------------------------------------'
    write(6,*) 'Monoexcitations to consider are: Excitation, Occupied Orbital, Virtual Orbital'

    call GenerateMonoexcitations

    do l1=1,nmonoex
        write(6,'(3I6)') l1,monoexcitations(l1,1),monoexcitations(l1,2)
    enddo

    if (uhf) then
        stop 'Stability analysis not implemented for UHF solutions'
    endif

    ! construct a and b matrices for singlet and triplet instablities
    allocate(asinglet(nmonoex,nmonoex),stat=ierr)
    if (ierr.ne.0) then
        stop 'Error allocating asinglet'
    endif
    allocate(bsinglet(nmonoex,nmonoex),stat=ierr)
    if (ierr.ne.0) then
        stop 'Error allocating bsinglet'
    endif
    allocate(atriplet(nmonoex,nmonoex),stat=ierr)
    if (ierr.ne.0) then
        stop 'Error allocating atriplet'
    endif
    allocate(btriplet(nmonoex,nmonoex),stat=ierr)
    if (ierr.ne.0) then
        stop 'Error allocating btriplet'
    endif

    write(6,*) '---------------------------------------------------------------'
    write(6,*) 'Constructing matrices A^1,A^3,B^1,B^3...'

    ! atriplet
    atriplet(:,:) = 0.0_dp
    do l1=1,nmonoex
        do l2=1,nmonoex
            atriplet(l1,l2) = -umat(monoexcitations(l1,2),monoexcitations(l2,1),&
                &monoexcitations(l2,2),monoexcitations(l1,1))
            if ((monoexcitations(l1,1).eq.monoexcitations(l2,1)).and.&
                &(monoexcitations(l1,2).eq.monoexcitations(l2,2))) then
                atriplet(l1,l2) = atriplet(l1,l2) + (arr(monoexcitations(l2,2)) -&
                    &arr(monoexcitations(l1,1)))
            endif
        enddo
    enddo

    ! btriplet
    btriplet(:,:) = 0.0_dp
    do l1=1,nmonoex
        do l2=1,nmonoex
            btriplet(l1,l2) = -umat(monoexcitations(l1,2),monoexcitations(l2,2),&
                &monoexcitations(l2,1),monoexcitations(l1,1))
        enddo
    enddo

    ! asinglet
    asinglet(:,:) = 0.0_dp
    do l1=1,nmonoex
        do l2=1,nmonoex
            asinglet(l1,l2) = atriplet(l1,l2) + (2.0_dp*&
                &umat(monoexcitations(l1,2),monoexcitations(l2,1),&
                &monoexcitations(l1,1),monoexcitations(l2,2)))
        enddo
    enddo

    ! bsinglet
    bsinglet(:,:) = 0.0_dp
    do l1=1,nmonoex
        do l2=1,nmonoex
            bsinglet(l1,l2) = btriplet(l1,l2) + (2.0_dp*&
                &umat(monoexcitations(l1,2),monoexcitations(l2,2),&
                &monoexcitations(l1,1),monoexcitations(l2,1)))
        enddo
    enddo

    ! since all orbitals are real, these can be combined to matrices
    allocate(xsinglet(nmonoex,nmonoex),stat=ierr)
    if (ierr.ne.0) then
        stop 'Error allocating xsinglet'
    endif
    allocate(xsinglet_breaktimerev(nmonoex,nmonoex),stat=ierr)
    if (ierr.ne.0) then
        stop 'Error allocating xsinglet_breaktimerev'
    endif
    allocate(xtriplet(nmonoex,nmonoex),stat=ierr)
    if (ierr.ne.0) then
        stop 'Error allocating xtriplet'
    endif
    allocate(xtriplet_breaktimerev(nmonoex,nmonoex),stat=ierr)
    if (ierr.ne.0) then
        stop 'Error allocating xtriplet_breaktimerev'
    endif
    allocate(xsinglet_eval(nmonoex),stat=ierr)
    if (ierr.ne.0) then
        stop 'Error allocating xsinglet'
    endif
    allocate(xsinglet_breaktimerev_eval(nmonoex),stat=ierr)
    if (ierr.ne.0) then
        stop 'Error allocating xsinglet_breaktimerev'
    endif
    allocate(xtriplet_eval(nmonoex),stat=ierr)
    if (ierr.ne.0) then
        stop 'Error allocating xtriplet'
    endif
    allocate(xtriplet_breaktimerev_eval(nmonoex),stat=ierr)
    if (ierr.ne.0) then
        stop 'Error allocating xtriplet_breaktimerev'
    endif
    lscr = 4*nmonoex
    allocate(scr(lscr),stat=ierr)
    if (ierr.ne.0) then
        stop 'Error allocating scr'
    endif
    allocate(inst_all(40,8),stat=ierr)
    if (ierr.ne.0) then
        stop 'Error allocating inst_xsinglet'
    endif
 
    inst_all(:,:) = 0

    ! Singlet instability which retains time reversal symmetry
    ! A1 + B1
    xsinglet(:,:) = 0.0_dp
    xsinglet = asinglet + bsinglet

    ! Triplet instability which retains time reversal symmetry
    ! A3 + B3
    xtriplet(:,:) = 0.0_dp
    xtriplet = atriplet + btriplet

    ! Singlet instability which breaks time reversal symmetry
    ! A1 - B1
    xsinglet_breaktimerev(:,:) = 0.0_dp
    xsinglet_breaktimerev = asinglet - bsinglet

    ! Triplet instability which breaks time reversal symmetry
    ! A3 - B3
    xtriplet_breaktimerev(:,:) = 0.0_dp
    xtriplet_breaktimerev = atriplet - btriplet
    

    ! Diagonalise instability matrice
    ! all matrices are hermitian
    ! since all orbitals are real, all matrices are symmetric

    write(6,*) 'Diagonalising xsinglet...'
    call dsyev('V','U',nmonoex,xsinglet,nmonoex,xsinglet_eval,scr,lscr,ierr)
    if (ierr.ne.0) then
        stop 'Error diagonalising xsinglet'
    endif
    write(6,*) 'Diagonalising xtriplet...'
    call dsyev('V','U',nmonoex,xtriplet,nmonoex,xtriplet_eval,scr,lscr,ierr)
    if (ierr.ne.0) then
        stop 'Error diagonalising xtinglet'
    endif
    write(6,*) 'Diagonalising xsinglet_breaktimerev...'
    call dsyev('V','U',nmonoex,xsinglet_breaktimerev,nmonoex,&
        &xsinglet_breaktimerev_eval,scr,lscr,ierr)
    if (ierr.ne.0) then
        stop 'Error diagonalising xsinglet_breaktimerev'
    endif
    write(6,*) 'Diagonalising xtriplet_breaktimerev...'
    call dsyev('V','U',nmonoex,xtriplet_breaktimerev,nmonoex,&
        &xtriplet_breaktimerev_eval,scr,lscr,ierr)
    if (ierr.ne.0) then
        stop 'Error diagonalising xtriplet_breaktimerev'
    endif

    write(6,*) '-------------------------------------------------------------'
    write(6,*) 'Eigenvalues of xsinglet'
    do l1=1,nmonoex
        write(6,*) l1,xsinglet_eval(l1)
    enddo
    write(6,*) '-------------------------------------------------------------'
    write(6,*) 'Eigenvalues of xtriplet'
    do l1=1,nmonoex
        write(6,*) l1,xtriplet_eval(l1)
    enddo
    write(6,*) '-------------------------------------------------------------'
    write(6,*) 'Eigenvalues of xsinglet_breaktimerev'
    do l1=1,nmonoex
        write(6,*) l1,xsinglet_breaktimerev_eval(l1)
    enddo
    write(6,*) '-------------------------------------------------------------'
    write(6,*) 'Eigenvalues of xtriplet_breaktimerev'
    do l1=1,nmonoex
        write(6,*) l1,xtriplet_breaktimerev_eval(l1)
    enddo

    ! the necessary and sufficient condition for the instability of the time reversal
    ! invariant closed-shell HF state is that either of these matrices
    ! has at least one negative eigenvalue

    write(6,*) '---------------------------------------------------------------'
    write(6,*) 'Searching for instabilities in HF solutions ...'
    write(6,*) '---------------------------------------------------------------'
    ! xsinglet
    m = 0
    do l1=1,nmonoex
        if (xsinglet_eval(l1).le.0.0_dp) then
            write(6,*) '--------------------------------------------------------'
            write(6,*) 'There is a singlet instability retaining time reversal &
                &symmetry with eigenvalue and eigenvector given by:'
            write(6,'(I7,5X,G20.12)') l1,xsinglet_eval(l1)
            do l2=1,nmonoex
                write(6,'(I5,G20.12)') l2,xsinglet(l2,l1)
            enddo
            m = m + 1
            inst_all(m,1) = l1
            inst_all(m,2) = xsinglet_eval(l1)
        endif
    enddo
    ! xtriplet
    m = 0
    do l1=1,nmonoex
        if (xtriplet_eval(l1).le.0.0_dp) then
            write(6,*) '-------------------------------------------------------'
            write(6,*) 'There is a triplet instability retaining time reversal &
                &symmetry with eigenvalue and eigenvector given by:'
            write(6,'(I7,5X,G20.12)') l1,xtriplet_eval(l1)
            do l2=1,nmonoex
                write(6,'(I5,G20.12)') l2,xtriplet(l2,l1)
            enddo
            m = m + 1
            inst_all(m,3) = l1
            inst_all(m,4) = xtriplet_eval(l1)
        endif
    enddo
    ! xsinglet_breaktimerev
    m = 0
    do l1=1,nmonoex
        if (xsinglet_breaktimerev_eval(l1).le.0.0_dp) then
            write(6,*) '----------------------------------------------------------'
            write(6,*) 'There is a singlet instability breaking time reversal &
                &symmetry with eigenvalue and eigenvector given by:'
            write(6,'(I7,5X,G20.12)') l1,xsinglet_breaktimerev_eval(l1)
            do l2=1,nmonoex
                write(6,'(I5,G20.12)') l2,xsinglet_breaktimerev(l2,l1)
            enddo
            m = m + 1
            inst_all(m,5) = l1
            inst_all(m,6) = xsinglet_breaktimerev_eval(l1)
        endif
    enddo
    ! xtriplet_breaktimerev
    m = 0
    do l1=1,nmonoex
        if (xtriplet_breaktimerev_eval(l1).le.0.0_dp) then
            write(6,*) '-------------------------------------------------------'
            write(6,*) 'There is a triplet instability breaking time reversal &
                &symmetry with eigenvalue and eigenvector given by:'
            write(6,'(I7,5X,G20.12)') l1,xtriplet_breaktimerev_eval(l1)
            do l2=1,nmonoex
                write(6,'(I5,G20.12)') l2,xtriplet_breaktimerev(l2,l1)
            enddo
            m = m + 1
            inst_all(m,7) = l1
            inst_all(m,8) = xtriplet_breaktimerev_eval(l1)
        endif
    enddo
 

    deallocate(asinglet)
    deallocate(atriplet)
    deallocate(bsinglet)
    deallocate(btriplet)

    ! find the type of instability which has the lowest eigenvalue
    ! dysev returns the eigenvalues in ascending order, so only
    ! the lowest one needs to be examined although
    ! others may be interesting as well ?

    allocate(temp_inds(norb,norb),stat=ierr)
    if (ierr.ne.0) then
        stop 'Error allocating temporary Umat'
    endif
    allocate(temp_tmat(norb,norb),stat=ierr)
    if (ierr.ne.0) then
        stop 'Error allocating temporary tmat'
    endif
    allocate(transform(norb,norb),stat=ierr)
    if (ierr.ne.0) then
        stop 'Error allocating transformation matrix'
    endif
    allocate(transarray(nmonoex),stat=ierr)
    if (ierr.ne.0) then
        stop 'Error allocating transarray'
    endif
 
    transform(:,:) = 0.0_dp
    temp_inds(:,:) = 0.0_dp
    temp_tmat(:,:) = 0.0_dp

    ! the orbitals of the symmetry-adapted RHF solution can be rotated
    ! so that the broken-symmetry solution is approximated/generated
    ! the ratios in which the virtual orbitals have to be admixed to the
    ! occupied HF orbitals in order to yield a variational wavefunction
    ! in the direaction of the mean-energy steepest descent are determined
    ! by the eigenvector of the stability analysis.
    ! a correct treatment involves an SCF procedure with the correct
    ! starting configuration
    ! to get a first approximation: it is assumed that each occupied
    ! orbital mixes together with only one virtual orbita
    ! this maintains the orthonormality using 

    ! A1+B1: x = a r + b r'
    ! A1-B1: x = a r + ib r' (a=a,b=-b)
    ! A3+B3: x = a r + b r' for a spins
    !        x = a r - n r' for b spins
    ! A3-B3: x = a r + ib r' for a spins
    !        x = a r - ib r' for b spins
    ! with a = cos a and b = sin a

    ! since only the rato is known one can calulate
    ! tan alpha 
    ! select lowest eigenvalue of all and 
    ! using tan b

    write(6,*) '------------------------------------------------------------'
    
    transarray(:) = 0.0_dp
    loc = minloc(inst_all)
    if (loc(2).eq.2) then
        write(6,*) 'Using singlet instability retaining time reversal symmetry'
        transarray(1:nmonoex) = xsinglet(1:nmonoex,int(inst_all(1,1)))
    elseif (loc(2).eq.4) then
        write(6,*) 'Using triplet instability retaining time reversal symmetry'
        transarray(1:nmonoex) = xtriplet(1:nmonoex,int(inst_all(1,3)))
    elseif (loc(2).eq.6) then
        write(6,*) 'Using singlet instability breaking time reversal symmetry'
        transarray(1:nmonoex) = xsinglet_breaktimerev(1:nmonoex,int(inst_all(1,5)))
    elseif (loc(2).eq.8) then
        write(6,*) 'Using triplet instability breaking time reversal symmetry'
        transarray(1:nmonoex) = xtriplet_breaktimerev(1:nmonoex,int(inst_all(1,7)))
    else
        write(6,*) 'There is no rotation to be performed'
    endif


    ! set up transformation matrix

    if (scheme.eq.1) then
        write(6,*) '---------------------------------------------------------------'
        write(6,*) 'Setting up transformation of orbitals via pairwise rotations...'
        !transform(:,:) = 0.0_dp
        !do l1=1,norb
        !    transform(l1,l1) = 1.0_dp
        !enddo
        do l1=1,nmonoex
            if (abs(transarray(l1)).gt.1e-8_dp) then
                transform(:,:) = 0.0_dp
                do l2=1,norb
                    transform(l2,l2) = 1.0_dp
                enddo
                ratio = transarray(l1)!/(1.0_dp-transarray(l1))
                angle = atan(ratio)
                i = monoexcitations(l1,1)
                a = monoexcitations(l1,2)
                transform(i,i) = cos(angle)
                transform(a,i) = sin(angle)
                transform(i,a) = -sin(angle)
                transform(a,a) = cos(angle)
                write(6,*) 'Rotations:',i,a,cos(angle),sin(angle)

                ! transform integrails
                ! dgemm performs the operation
                ! C = alpha * op(A) * op(B) + beta * op(c)
                ! alpha,beta: scalars
                ! op(a): m by k matrix
                ! op(b): k by n matrix
                ! op(C): m by n matrix
                do i=1,norb
                    do j=1,norb
                        ! (c_i\alpha)^T <ij|kl> = <\alpha j|kl>
                        temp_inds(:,:) = 0.0_dp
                        call dgemm('T','N',norb,norb,norb,1.0_dp,transform(:,:),norb,&
                        &umat(1:norb,1:norb,i,j),norb,0.0_dp,&
                        &temp_inds(1:norb,1:norb),norb)
                        ! (c_j\beta)^T (<\alpha j|kl>)^T
                        ! the transpose is needed
                        umat(1:norb,1:norb,i,j) = 0.0_dp
                        call dgemm('T','T',norb,norb,norb,1.0_dp,transform(:,:),norb,&
                        &temp_inds(1:norb,1:norb),norb,0.0_dp,&
                        &umat(1:norb,1:norb,i,j),norb)
                    enddo
                enddo
                do i=1,norb
                    do j=1,norb
                        ! (c_k\gamma)^T <\alpha \beta|kl> = <\alpha \beta|\gammal>
                        temp_inds(:,:) = 0.0_dp
                        call dgemm('T','N',norb,norb,norb,1.0_dp,transform(:,:),norb,&
                        &umat(i,j,1:norb,1:norb),norb,0.0_dp,&
                        &temp_inds(1:norb,1:norb),norb)
                        ! (c_l\delta)^T (<\alpha \beta|\gamma l>)^T = <\alpha \beta|\gamma l>
                        ! the transpose is needed
                        umat(i,j,1:norb,1:norb) = 0.0_dp
                        call dgemm('T','T',norb,norb,norb,1.0_dp,transform(:,:),norb,&
                        &temp_inds(1:norb,1:norb),norb,0.0_dp,&
                        &umat(i,j,1:norb,1:norb),norb)
                    enddo
                enddo

         
                ! since the transform one electron integrals
                ! tmat needs to be transformed as well
                temp_tmat(:,:) = 0.0_dp
                ! h_ij c_j\alpha = h_i\alpha
                call dgemm('N','N',norb,norb,norb,1.0_dp,tmat,norb,&
                &transform,norb,0.0_dp,temp_tmat,norb)
                ! (c_i\beta)^T h_i\alpha = h_\beta\alpha
                tmat(:,:) = 0.0_dp
                call dgemm('T','N',norb,norb,norb,1.0_dp,transform,norb,&
                &temp_tmat,norb,0.0_dp,tmat,norb)

            endif
        enddo

        !orbital energies
        call DiagFockElements(fdiag)
        arr = fdiag

        ! new HF ground state energy approximately since
        ! fock matrix is no longer diagonal

        hf_energy = 0.0_dp
        if (uhf) then
            do l1=1,nelec
                hf_energy = hf_energy + arr(l1)
            enddo
        elseif (.not.uhf) then
            do l1=1,(nelec/2)
                hf_energy = hf_energy + arr(l1)
            enddo
        endif

        write(6,*) '---------------------------------------------------'
        write(6,*) 'The new HF ground state energy is:',hf_energy
        write(6,*) '---------------------------------------------------'

    elseif (scheme.eq.2) then
        ! Perform a transformation using the vector connected to the
        ! lowest eigenvalue
        write(6,*) '---------------------------------------------------------------------'
        write(6,*) 'Setting up transformation of orbitals via eigenvector and symmetric Loewdin orthogonalisation...'
        transform(:,:) = 0.0_dp
        do l1=1,norb
            transform(l1,l1) = 1.0_dp
        enddo
        do l1=1,nmonoex
            if (abs(transarray(l1)).gt.1e-8_dp) then
                ratio = transarray(l1)!/(1.0_dp-transarray(l1))
                i = monoexcitations(l1,1)
                a = monoexcitations(l1,2)
                transform(i,i) = 1.0_dp
                transform(a,i) = ratio
                transform(i,a) = -ratio
                transform(a,a) = 1.0_dp
                write(6,*) 'Rotation:',i,a,transform(a,i),transform(i,a)
            endif
        enddo

        ! transform integrails
        ! dgemm performs the operation
        ! C = alpha * op(A) * op(B) + beta * op(c)
        ! alpha,beta: scalars
        ! op(a): m by k matrix
        ! op(b): k by n matrix
        ! op(C): m by n matrix
        do i=1,norb
            do j=1,norb
                ! (c_i\alpha)^T <ij|kl> = <\alpha j|kl>
                temp_inds(:,:) = 0.0_dp
                call dgemm('T','N',norb,norb,norb,1.0_dp,transform(:,:),norb,&
                &umat(1:norb,1:norb,i,j),norb,0.0_dp,&
                &temp_inds(1:norb,1:norb),norb)
                ! (c_j\beta)^T (<\alpha j|kl>)^T
                ! the transpose is needed
                umat(1:norb,1:norb,i,j) = 0.0_dp
                call dgemm('T','T',norb,norb,norb,1.0_dp,transform(:,:),norb,&
                &temp_inds(1:norb,1:norb),norb,0.0_dp,&
                &umat(1:norb,1:norb,i,j),norb)
            enddo
        enddo
        do i=1,norb
            do j=1,norb
                ! (c_k\gamma)^T <\alpha \beta|kl> = <\alpha \beta|\gammal>
                temp_inds(:,:) = 0.0_dp
                call dgemm('T','N',norb,norb,norb,1.0_dp,transform(:,:),norb,&
                &umat(i,j,1:norb,1:norb),norb,0.0_dp,&
                &temp_inds(1:norb,1:norb),norb)
                ! (c_l\delta)^T (<\alpha \beta|\gamma l>)^T = <\alpha \beta|\gamma l>
                ! the transpose is needed
                umat(i,j,1:norb,1:norb) = 0.0_dp
                call dgemm('T','T',norb,norb,norb,1.0_dp,transform(:,:),norb,&
                &temp_inds(1:norb,1:norb),norb,0.0_dp,&
                &umat(i,j,1:norb,1:norb),norb)
            enddo
        enddo

 
        ! since the transform one electron integrals
        ! tmat needs to be transformed as well
        temp_tmat(:,:) = 0.0_dp
        ! h_ij c_j\alpha = h_i\alpha
        call dgemm('N','N',norb,norb,norb,1.0_dp,tmat,norb,&
        &transform,norb,0.0_dp,temp_tmat,norb)
        ! (c_i\beta)^T h_i\alpha = h_\beta\alpha
        tmat(:,:) = 0.0_dp
        call dgemm('T','N',norb,norb,norb,1.0_dp,transform,norb,&
        &temp_tmat,norb,0.0_dp,tmat,norb)


        ! perform a Loewdin symmetric orthogonalisation
        ! this gives the basis which resembles most closely the
        ! original orbitals

        ! overlap matrix
        allocate(smat(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating smat'
        endif
        allocate(smat_eval(norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating smat_eval'
        endif
        allocate(temp_smat(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating temp_smat'
        endif
        allocate(temp2_smat(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating temp_smat'
        endif
 
        smat(:,:) = 0.0_dp
        do l1=1,norb
            do l2=1,norb
                ! the hf orbitals are orthonormal
                ! hence S_ij = sum_{k} (c_ki)*(c_kj)
                do l3=1,norb
                    smat(l2,l1) = smat(l2,l1) + (transform(l3,l2)*&
                        &transform(l3,l1))
                enddo
            enddo
        enddo

        ! using S^-1/2 = A sqrt(B) A^T
        ! where S = A B A^T
        ! and B is the diagonal matrix with e^-1/2
        ! where e are the eigenvalues of S
        call dsyev('V','U',norb,smat,norb,smat_eval,scr,lscr,ierr)
        if (ierr.ne.0) then
            stop 'Error diagonalising smat'
        endif

        temp_smat(:,:) = 0.0_dp
        do l1=1,norb
            temp_smat(l1,l1) = 1.0_dp/sqrt(smat_eval(l1))
        enddo
        temp2_smat(:,:) = 0.0_dp
        ! temp_smat_ij (smat_kj)^T = temp2_smat_ik
        call dgemm('N','T',norb,norb,norb,1.0_dp,temp_smat,norb,&
        &smat(1:norb,1:norb),norb,0.0_dp,&
        &temp2_smat(1:norb,1:norb),norb)
        ! smat_ij temp2_smat_jk = temp_smat_ik 
        temp_smat(:,:) = 0.0_dp
        call dgemm('N','N',norb,norb,norb,1.0_dp,smat(:,:),norb,&
        &temp2_smat(1:norb,1:norb),norb,0.0_dp,&
        &temp_smat(1:norb,1:norb),norb)


        ! transform integrails
        ! dgemm performs the operation
        ! C = alpha * op(A) * op(B) + beta * op(c)
        ! alpha,beta: scalars
        ! op(a): m by k matrix
        ! op(b): k by n matrix
        ! op(C): m by n matrix
        do i=1,norb
            do j=1,norb
                ! (c_i\alpha)^T <ij|kl> = <\alpha j|kl>
                temp_inds(:,:) = 0.0_dp
                call dgemm('T','N',norb,norb,norb,1.0_dp,temp_smat(:,:),norb,&
                &umat(1:norb,1:norb,i,j),norb,0.0_dp,&
                &temp_inds(1:norb,1:norb),norb)
                ! (c_j\beta)^T (<\alpha j|kl>)^T
                ! the transpose is needed
                umat(1:norb,1:norb,i,j) = 0.0_dp
                call dgemm('T','T',norb,norb,norb,1.0_dp,temp_smat(:,:),norb,&
                &temp_inds(1:norb,1:norb),norb,0.0_dp,&
                &umat(1:norb,1:norb,i,j),norb)
            enddo
        enddo
        do i=1,norb
            do j=1,norb
                ! (c_k\gamma)^T <\alpha \beta|kl> = <\alpha \beta|\gammal>
                temp_inds(:,:) = 0.0_dp
                call dgemm('T','N',norb,norb,norb,1.0_dp,temp_smat(:,:),norb,&
                &umat(i,j,1:norb,1:norb),norb,0.0_dp,&
                &temp_inds(1:norb,1:norb),norb)
                ! (c_l\delta)^T (<\alpha \beta|\gamma l>)^T = <\alpha \beta|\gamma l>
                ! the transpose is needed
                umat(i,j,1:norb,1:norb) = 0.0_dp
                call dgemm('T','T',norb,norb,norb,1.0_dp,temp_smat(:,:),norb,&
                &temp_inds(1:norb,1:norb),norb,0.0_dp,&
                &umat(i,j,1:norb,1:norb),norb)
            enddo
        enddo

 
        ! since the transform one electron integrals
        ! tmat needs to be transformed as well
        temp_tmat(:,:) = 0.0_dp
        ! h_ij c_j\alpha = h_i\alpha
        call dgemm('N','N',norb,norb,norb,1.0_dp,temp_smat,norb,&
        &transform,norb,0.0_dp,temp_tmat,norb)
        ! (c_i\beta)^T h_i\alpha = h_\beta\alpha
        tmat(:,:) = 0.0_dp
        call dgemm('T','N',norb,norb,norb,1.0_dp,temp_smat,norb,&
        &temp_tmat,norb,0.0_dp,tmat,norb)



        !orbital energies
        call DiagFockElements(fdiag)
        arr = fdiag

        ! new HF ground state energy approximately since
        ! fock matrix is no longer diagonal

        hf_energy = 0.0_dp
        if (uhf) then
            do l1=1,nelec
                hf_energy = hf_energy + arr(l1)
            enddo
        elseif (.not.uhf) then
            do l1=1,(nelec/2)
                hf_energy = hf_energy + arr(l1)
            enddo
        endif

        write(6,*) '---------------------------------------------------'
        write(6,*) 'The new HF ground state energy is:',hf_energy
        write(6,*) '---------------------------------------------------'

        deallocate(temp_smat)
        deallocate(temp2_smat)
        deallocate(smat)
 
    elseif (scheme.eq.3) then

        ! perform RHF on starting guess wavefunction
        ! need complex orbitals
        write(6,*) '---------------------------------------------------------------------'
        write(6,*) 'Setting up transformation of orbitals via diagonalising Fock matrix...'
        allocate(coeffs(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating coeffs'
        endif
        allocate(grade1(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating grade1'
        endif
        allocate(temptrans(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating temptrans'
        endif
        allocate(grade2(norb,norb,norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating grade2'
        endif

        ! generate starting wavefunction using 'reverse perturbation'

        !call StartingWavefunction(coeffs)

        ! the calculated instability matrices correspond to the Hessian matrices
        ! the eigenvector corresponding to the smallest eigenvalue will give
        ! the gradient in the direction of steepest descent

        write(6,*) 'The gradient obtained from the instability analysis is: E^(1)_pq'
        transform(:,:) = 0.0_dp
        do l1=1,norb
            transform(l1,l1) = 1.0_dp
        enddo
        do l1=1,nmonoex
            if (abs(transarray(l1)).gt.1e-8_dp) then
                ratio = transarray(l1)!/(1.0_dp-transarray(l1))
                i = monoexcitations(l1,1)
                a = monoexcitations(l1,2)
                transform(i,i) = 1.0_dp
                transform(a,i) = ratio
                transform(i,a) = -ratio
                transform(a,a) = 1.0_dp
            endif
        enddo

        do l1=1,norb
            do l2=1,norb
                if (abs(transform(l2,l1)).gt.1e-12_dp) then
                    write(6,*) l2,l1,transform(l2,l1)
                endif
            enddo
        enddo

        coeffs = cmplx(transform,0.0_dp,dp)
 
        ! this wourine generates a starting guess wavefunction using a
        ! 'reverse perturbation' approach 
        !call StartingWavefunction(coeffs)

        ! This combination of routines would perform a steepest descent type
        ! implementation followed by an RHF procedure
        ! follow a steepest descent way d
        !call SteepestDescent(coeffs,grade1,grade2)

        ! these two routines would perform a standart RHF procedure either
        ! without or with DIIS
        ! the DIIS approach tends to lead back to the original solution
        !call RHF(coeffs,fdiag) 
        !call RHF_DIIS(coeffs,fdiag)

        
        ! this uses a Newton-Raphson iteration scheme in order to
        ! find a new minimum based on the previous instability analysis
        grade1 = coeffs
        temptrans(:,:) = 0.0_dp
        call TransformFromGradient(temptrans,grade1,1.0_dp)
 
        grade1(:,:) = 0.0_dp
        do l1=1,norb
            grade1(l1,l1) = 1.0_dp
        enddo

        coeffs(:,:) = 0.0_dp
        call zgemm('T','N',norb,norb,norb,(1.0_dp,0.0_dp),temptrans,norb,&
        &grade1,norb,(0.0_dp,0.0_dp),coeffs,norb)

        !grade1(:,:) = 0.0_dp
        grade2(:,:,:,:) = 0.0_dp

        call PerformNewtonRaphson(coeffs,grade1,grade2)

        ! this diagonalises the fock matrix just once
        !call RHF(coeffs,fdiag)


        ! as long as all coefficients are real this will be okay
        transform = real(coeffs,dp)
        arr = fdiag

        deallocate(coeffs)
        deallocate(temptrans)
        deallocate(grade1)
        deallocate(grade2)

        ! transform integrails
        ! dgemm performs the operation
        ! C = alpha * op(A) * op(B) + beta * op(c)
        ! alpha,beta: scalars
        ! op(a): m by k matrix
        ! op(b): k by n matrix
        ! op(C): m by n matrix
        do i=1,norb
            do j=1,norb
                ! (c_i\alpha)^T <ij|kl> = <\alpha j|kl>
                temp_inds(:,:) = 0.0_dp
                call dgemm('T','N',norb,norb,norb,1.0_dp,transform(:,:),norb,&
                &umat(1:norb,1:norb,i,j),norb,0.0_dp,&
                &temp_inds(1:norb,1:norb),norb)
                ! (c_j\beta)^T (<\alpha j|kl>)^T
                ! the transpose is needed
                umat(1:norb,1:norb,i,j) = 0.0_dp
                call dgemm('T','T',norb,norb,norb,1.0_dp,transform(:,:),norb,&
                &temp_inds(1:norb,1:norb),norb,0.0_dp,&
                &umat(1:norb,1:norb,i,j),norb)
            enddo
        enddo
        do i=1,norb
            do j=1,norb
                ! (c_k\gamma)^T <\alpha \beta|kl> = <\alpha \beta|\gammal>
                temp_inds(:,:) = 0.0_dp
                call dgemm('T','N',norb,norb,norb,1.0_dp,transform(:,:),norb,&
                &umat(i,j,1:norb,1:norb),norb,0.0_dp,&
                &temp_inds(1:norb,1:norb),norb)
                ! (c_l\delta)^T (<\alpha \beta|\gamma l>)^T = <\alpha \beta|\gamma l>
                ! the transpose is needed
                umat(i,j,1:norb,1:norb) = 0.0_dp
                call dgemm('T','T',norb,norb,norb,1.0_dp,transform(:,:),norb,&
                &temp_inds(1:norb,1:norb),norb,0.0_dp,&
                &umat(i,j,1:norb,1:norb),norb)
            enddo
        enddo

        ! since the transform one electron integrals
        ! tmat needs to be transformed as well
        temp_tmat(:,:) = 0.0_dp
        ! h_ij c_j\alpha = h_i\alpha
        call dgemm('N','N',norb,norb,norb,1.0_dp,tmat,norb,&
        &transform,norb,0.0_dp,temp_tmat,norb)
        ! (c_i\beta)^T h_i\alpha = h_\beta\alpha
        tmat(:,:) = 0.0_dp
        call dgemm('T','N',norb,norb,norb,1.0_dp,transform,norb,&
        &temp_tmat,norb,0.0_dp,tmat,norb)



        ! new HF ground state energy approximately since
        ! fock matrix is no longer diagonal

        hf_energy = 0.0_dp
        if (uhf) then
            do l1=1,nelec
                hf_energy = hf_energy + arr(l1) + tmat(l1,l1)
            enddo
        elseif (.not.uhf) then
            do l1=1,(nelec/2)
                hf_energy = hf_energy + arr(l1) + tmat(l1,l1)
            enddo
        endif

        write(6,*) '---------------------------------------------------'
        write(6,*) 'The new HF ground state energy is:',hf_energy
        write(6,*) '---------------------------------------------------'


    endif


    deallocate(xsinglet)
    deallocate(xsinglet_eval)
    deallocate(xtriplet)
    deallocate(xtriplet_eval)
    deallocate(xsinglet_breaktimerev)
    deallocate(xsinglet_breaktimerev_eval)
    deallocate(xtriplet_breaktimerev)
    deallocate(xtriplet_breaktimerev_eval)
    deallocate(inst_all)
    deallocate(temp_inds)
    deallocate(temp_tmat)
    deallocate(scr)

    ! writing out a new FCIUDMP_file file with rotated orbital integrals

    write(6,*) 'Writing FCIDUMP_rotated file...'

    call WriteFCIDUMP

    deallocate(umat)
    deallocate(tmat)
    deallocate(arr)
    deallocate(transform)
    deallocate(transarray)
    deallocate(monoexcitations)
 
contains


    subroutine WriteFCIDUMP
 
        ! this subroutine writes out a new FCIUDMP_file file 
        ! with rotated orbital integrals
        ! the orbital ordering and format is the same as in
        ! the original FCIDUMP file


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


    end subroutine WriteFCIDUMP


    subroutine GenerateMonoexcitations

        ! this subroutine is to generate all possible monoexcitations from
        ! the occupied orbitals to the virtual orbitals


        integer :: l1,l2,m
        integer :: g,b


        if (tMolpro) then
            if (uhf) then
                ! spin orbitals
                m = 0
                do l1=1,nelec
                    do l2=(nelec+1),norb
                        m = m + 1
                        ! occupied orbital
                        monoexcitations(m,1) = l1
                        ! virtual orbital
                        monoexcitations(m,2) = l2
                    enddo
                enddo
            elseif (.not.uhf) then
                ! spatial orbitals
                m = 0
                do l1=1,(nelec/2)
                    do l2=((nelec/2)+1),norb
                        m = m + 1
                        ! occupied orbital
                        monoexcitations(m,1) = l1
                        ! virtual orbital
                        monoexcitations(m,2) = l2
                    enddo
                enddo
            endif
        else
            if (uhf) then
                ! spin orbitals
                m = 0
                do l1=1,nelec
                    g = energyorder(l1)
                    do l2=(nelec+1),norb
                        b = energyorder(l2)
                        m = m + 1
                        ! occupied orbital
                        monoexcitations(m,1) = g
                        ! virtual orbital
                        monoexcitations(m,2) = b
                    enddo
                enddo
            elseif (.not.uhf) then
                ! spatial orbitals
                m = 0
                do l1=1,(nelec/2)
                    g = energyorder(l1)
                    do l2=((nelec/2)+1),norb
                        m = m + 1
                        b = energyorder(l2)
                        ! occupied orbital
                        monoexcitations(m,1) = g
                        ! virtual orbital
                        monoexcitations(m,2) = b
                    enddo
                enddo
            endif
        endif

        if (m.ne.nmonoex) then
            stop 'Error: monoexcitations are not enumerated correctly'
        endif



    end subroutine GenerateMonoexcitations


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


    subroutine RHF(coeffs,evals)

        ! This subroutine performs an RHF calculation on a given starting set
        ! of coefficients in order to find a different solution

        complex(dp), allocatable, intent(inout) :: coeffs(:,:)
        real(dp), allocatable, intent(inout) :: evals(:)
        complex(dp), allocatable :: fmat(:,:),dmat(:,:),dmat_old(:,:)
        real(dp), allocatable :: rrscr(:),rscr(:)
        complex(dp) :: energy_old,energy,std
        logical :: ints_complex
        integer :: rlscr,l1,l2,l3,l4,ierr,iter
        
        rlscr = 4*norb
        allocate(fmat(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating fmat'
        endif
        allocate(dmat(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating dmat'
        endif
        allocate(dmat_old(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating dmat_old'
        endif
        allocate(rrscr(rlscr),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating rrscr'
        endif
        allocate(rscr(rlscr),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating rscr'
        endif

        write(6,*) '--------------------------------------------------'
        write(6,*) 'Performing RHF SCF on given starting wavefunction'
        write(6,*) '--------------------------------------------------'

        ! construct density matrix
        dmat(:,:) = 0.0_dp
        do l1=1,norb
            do l2=1,norb
                do l3=1,(nelec/2)
                    dmat(l2,l1) = dmat(l2,l1) + (2.0_dp*coeffs(l2,l3)*&
                        &conjg(coeffs(l1,l2)))
                enddo
            enddo
        enddo

        ! fock matrix
        fmat(:,:) = 0.0_dp
        do l1=1,norb
            do l2=1,norb
                do l3=1,norb
                    do l4=1,norb
                        fmat(l2,l1) = fmat(l2,l1) + (dmat(l4,l3)*&
                            &(umat(l2,l3,l1,l4)-(0.5_dp*umat(l2,l3,l4,l1))))
                    enddo
                enddo
            enddo
        enddo

        fmat = fmat + tmat

        iter = 0
        energy = 0.0_dp
        energy_old = 0.0_dp
        std = 0.0_dp
        dmat_old = dmat

        do
            iter = iter + 1
            write(6,*) '----------------------------------------------'
            write(6,*) 'Iteration:',iter
            write(6,*) '----------------------------------------------'


    !        ! mix in part of previous density matrix to avoid oscillation
    !        ! between solution
    !        if (.true.) then!.and.(mod(iter,2).eq.0)) then

    !            dmat = (0.95_dp*dmat_old) + ((1.0_dp-0.95_dp)*dmat)
    !            ! fock matrix
    !            fmat(:,:) = 0.0_dp
    !            do l1=1,norb
    !                do l2=1,norb
    !                    do l3=1,norb
    !                        do l4=1,norb
    !                            fmat(l2,l1) = fmat(l2,l1) + (dmat(l4,l3)*&
    !                                &(umat(l2,l3,l1,l4)-(0.5_dp*umat(l2,l3,l4,l1))))
    !                        enddo
    !                    enddo
    !                enddo
    !            enddo

    !            fmat = fmat + tmat

    !        endif

            ! diagonalise fmat
            write(6,*) 'Diagonalising fmat...'
            call zheev('V','U',norb,fmat,norb,evals,rscr,rlscr,rrscr,ierr)
            if (ierr.ne.0) then
                stop 'Error diagonalising fmat'
            endif

            write(6,*) 'Eigenvectors:'
            do l1=1,norb
                write(6,'(i20)') l1
                do l2=1,norb
                    write(6,*) fmat(l2,l1)
                enddo
            enddo

            write(6,*) 'Eigenvalues'
            do l1=1,norb
                write(6,*) l1,evals(l1)
            enddo

            coeffs(:,:) = fmat(:,:)

            ! construct density matrix
            dmat_old = dmat
            dmat(:,:) = 0.0_dp
            do l1=1,norb
                do l2=1,norb
                    do l3=1,(nelec/2)
                        dmat(l2,l1) = dmat(l2,l1) + (2.0_dp*coeffs(l2,l3)*&
                            &conjg(coeffs(l1,l2)))
                    enddo
                enddo
            enddo

            ! fock matrix
            fmat(:,:) = 0.0_dp
            do l1=1,norb
                do l2=1,norb
                    do l3=1,norb
                        do l4=1,norb
                            fmat(l2,l1) = fmat(l2,l1) + (dmat(l4,l3)*&
                                &(umat(l2,l3,l1,l4)-(0.5_dp*umat(l2,l3,l4,l1))))
                        enddo
                    enddo
                enddo
            enddo

            fmat = fmat + tmat

    
            ! ground state energy and difference in density matrix as 
            ! convergence criterion
            energy_old = energy
            energy = 0.0_dp
            do l1=1,norb
                do l2=1,norb
                    energy = energy + (0.5_dp*dmat(l1,l2)*&
                        &(tmat(l2,l1)+fmat(l2,l1)))
                enddo
            enddo

            ! difference in density matrices
            std = 0.0_dp
            do l1=1,norb
                do l2=1,norb
                    std = std + (abs(dmat(l2,l1)-dmat_old(l2,l1)))
                enddo
            enddo

            if (abs(aimag(energy)).gt.1e-12_dp) then
                stop 'Error: energy is complex'
            endif

            write(6,*) '-------------------------------------------------'
            write(6,*) 'The HF ground state energy of the previous and &
                & current iterations are'
            write(6,*) energy_old,energy
            write(6,*) '--------------------------------------------------'
            write(6,*) 'The summed differences in the density matrices of the &
                &previous and current iterations are'
            write(6,*) std
            write(6,*) '---------------------------------------------------'

            ! only diagonalise fock matrix one in order to move one 
            ! step into the correct direction
            if (iter.eq.1) then
            !if ((abs(energy-energy_old).lt.1e-12_dp).and.&
            !    &(abs(std).lt.1e-12_dp)) then
                write(6,*) 'Convergence achieved: stopping calculation !'
                write(6,*) '---------------------------------------------'
                exit
            endif
        enddo


        ! test whether coefficients are complex or not
        ints_complex = .false.
        do l1=1,norb
            do l2=1,norb
                if (aimag(coeffs(l2,l1)).gt.1e-12_dp) then
                    ints_complex = .true.
                endif
            enddo
        enddo

        if (ints_complex) then
            write(6,*) 'Coefficients are complex'
        else
            write(6,*) 'Coefficients and hence integrals are all real'
        endif


        deallocate(fmat)
        deallocate(dmat)
        deallocate(dmat_old)
        deallocate(rscr)
        deallocate(rrscr)

    end subroutine RHF

    
    subroutine StartingWavefunction(coeffs)

        ! this subroutine is to find a good starting wavefunction for the RHF
        ! calculation in order to converge to the broken symmetry solution
        ! this employs a 'reverse perturbation' of the fock matrix using 
        ! the eigenvectors of the instablity problem

        complex(dp), allocatable, intent(inout) :: coeffs(:,:)
        real(dp), allocatable :: hfevals(:)
        complex(dp), allocatable :: fmat(:,:),dmat(:,:)
        real(dp), allocatable :: srrscr(:),srscr(:)
        complex(dp) :: energy
        real(dp), allocatable :: lambda(:),hf_energy(:)
        integer :: srlscr,l1,l2,l4,l5,ierr
        integer :: minlambda(1)
        
        srlscr = 4*norb
        allocate(fmat(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating fmat'
        endif
        allocate(dmat(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating dmat'
        endif
        allocate(srrscr(srlscr),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating srrscr'
        endif
        allocate(srscr(lscr),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating srscr'
        endif
        allocate(lambda(30),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating lambda'
        endif
        allocate(hf_energy(30),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating hf_energy'
        endif
        allocate(hfevals(norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating hfevals'
        endif



        write(6,*) '--------------------------------------------------'
        write(6,*) 'Finding starting wavefunction using "reverse perturbation"'
        write(6,*) '--------------------------------------------------'

        lambda(1) = 0.0_dp
        do l5=2,30
            lambda(l5) = lambda(l5-1) + (1.0_dp/30.0_dp)
        enddo

        do l5=1,30
            coeffs(:,:) = 1.0_dp
            
            ! construct density matrix
            dmat(:,:) = 0.0_dp
            do l1=1,norb
                do l2=1,norb
                    do l3=1,(nelec/2)
                        dmat(l2,l1) = dmat(l2,l1) + (2.0_dp*coeffs(l2,l3)*&
                            &conjg(coeffs(l1,l2)))
                    enddo
                enddo
            enddo

            ! fock matrix
            fmat(:,:) = 0.0_dp
            do l1=1,norb
                do l2=1,norb
                    do l3=1,norb
                        do l4=1,norb
                            fmat(l2,l1) = fmat(l2,l1) + (dmat(l4,l3)*&
                            &(umat(l2,l3,l1,l4)-(0.5_dp*umat(l2,l3,l4,l1))))
                        enddo
                    enddo
                enddo
            enddo

            fmat = fmat + tmat

            do l1=1,nmonoex
                ratio = transarray(l1)
                i = monoexcitations(l1,1)
                a = monoexcitations(l1,2)
                fmat(i,a) = fmat(i,a) + (lambda(l5)*ratio*(arr(a)-arr(i)))
            enddo
 
 
            ! diagonalise fmat
            !write(6,*) 'Diagonalising fmat...'
            call zheev('V','U',norb,fmat,norb,hfevals,srscr,srlscr,srrscr,ierr)
            if (ierr.ne.0) then
                stop 'Error diagonalising fmat'
            endif

        !    write(6,*) 'Eigenvectors:'
        !    do l1=1,norb
        !        write(6,'(i20)') l1
        !        do l2=1,norb
        !            write(6,*) fmat(l2,l1)
        !        enddo
        !    enddo

        !    write(6,*) 'Eigenvalues'
        !    do l1=1,norb
        !        write(6,*) l1,evals(l1)
        !    enddo

            coeffs(:,:) = fmat(:,:)
            fdiag(:) = hfevals(:)

            ! construct density matrix
            dmat(:,:) = 0.0_dp
            do l1=1,norb
                do l2=1,norb
                    do l3=1,(nelec/2)
                        dmat(l2,l1) = dmat(l2,l1) + (2.0_dp*coeffs(l2,l3)*&
                            &conjg(coeffs(l1,l2)))
                    enddo
                enddo
            enddo

            ! fock matrix
            fmat(:,:) = 0.0_dp
            do l1=1,norb
                do l2=1,norb
                    do l3=1,norb
                        do l4=1,norb
                            fmat(l2,l1) = fmat(l2,l1) + (dmat(l4,l3)*&
                                &(umat(l2,l3,l1,l4)-(0.5_dp*umat(l2,l3,l4,l1))))
                        enddo
                    enddo
                enddo
            enddo

            fmat = fmat + tmat

    
            ! ground state energy and difference in density matrix as 
            ! convergence criterion
            energy = 0.0_dp
            do l1=1,norb
                do l2=1,norb
                    energy = energy + (0.5_dp*dmat(l1,l2)*&
                        &(tmat(l2,l1)+fmat(l2,l1)))
                enddo
            enddo

            if (abs(aimag(energy)).gt.1e-12_dp) then
                stop 'Error: energy is complex'
            endif

            hf_energy(l5) = real(energy,dp)

        enddo

        write(6,*) '--------------------------------------------------------------------'
        write(6,*) 'The pertrubation parameter lambda and HF ground state energies are'
        do l5=1,30
            write(6,*) l5,lambda(l5),hf_energy(l5)
        enddo

        ! choose the lambe which gives the lowest hf_energy
        minlambda = minloc(hf_energy)
        write(6,*) '-------------------------------------------------------------------'
        write(6,*) 'Choosing:',lambda(minlambda(1)),hf_energy(minlambda(1))
        write(6,*) '-------------------------------------------------------------------'

        ! evaluate coefficients for a starting wavefunction
         coeffs(:,:) = 1.0_dp
        
        ! construct density matrix
        dmat(:,:) = 0.0_dp
        do l1=1,norb
            do l2=1,norb
                do l3=1,(nelec/2)
                    dmat(l2,l1) = dmat(l2,l1) + (2.0_dp*coeffs(l2,l3)*&
                        &conjg(coeffs(l1,l2)))
                enddo
            enddo
        enddo

        ! fock matrix
        fmat(:,:) = 0.0_dp
        do l1=1,norb
            do l2=1,norb
                do l3=1,norb
                    do l4=1,norb
                        fmat(l2,l1) = fmat(l2,l1) + (dmat(l4,l3)*&
                        &(umat(l2,l3,l1,l4)-(0.5_dp*umat(l2,l3,l4,l1))))
                    enddo
                enddo
            enddo
        enddo

        fmat = fmat + tmat

        do l1=1,nmonoex
            ratio = transarray(l1)
            i = monoexcitations(l1,1)
            a = monoexcitations(l1,2)
            fmat(i,a) = fmat(i,a) + (lambda(minlambda(1))*ratio*(arr(a)-arr(i)))
        enddo
 
 
        ! diagonalise fmat
        !write(6,*) 'Diagonalising fmat...'
        call zheev('V','U',norb,fmat,norb,hfevals,srscr,srlscr,srrscr,ierr)
        if (ierr.ne.0) then
            stop 'Error diagonalising fmat'
        endif

        coeffs = fmat

     !   do l1=1,norb
     !       do l2=1,norb
     !           write(6,*) l1,l2,coeffs(l2,l2)
     !       enddo
     !   enddo



    end subroutine StartingWavefunction

    subroutine RHF_DIIS(coeffs,evals)

        ! This subroutine performs an RHF calculation on a given starting set
        ! of coefficients in order to find a different solution

        complex(dp), allocatable, intent(inout) :: coeffs(:,:)
        real(dp), allocatable, intent(inout) :: evals(:)
        complex(dp), allocatable :: fmat(:,:),dmat(:,:),dmat_old(:,:)
        real(dp), allocatable :: rrscr(:),rscr(:)
        complex(dp) :: energy_old,energy,std
        logical :: ints_complex
        integer :: rlscr,l1,l2,l3,l4,ierr,iter,nerror,beforeiter
        complex(dp), allocatable :: f1(:,:),f2(:,:),f3(:,:),f4(:,:),f5(:,:)
        complex(dp), allocatable :: f6(:,:),f7(:,:),f8(:,:)
        complex(dp), allocatable :: e1(:,:),e2(:,:),e3(:,:),e4(:,:),e5(:,:)
        complex(dp), allocatable :: e6(:,:),e7(:,:),e8(:,:)
        complex(dp), allocatable :: bmat(:,:),bmat_coeffs(:)
        real(dp), allocatable :: bmat_evals(:)
        integer, allocatable :: ipiv(:)

        rlscr = 4*norb
        allocate(fmat(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating fmat'
        endif
        allocate(dmat(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating dmat'
        endif
        allocate(dmat_old(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating dmat_old'
        endif
        allocate(rrscr(rlscr),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating rrscr'
        endif
        allocate(rscr(rlscr),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating rscr'
        endif
        allocate(e1(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating e1'
        endif
        allocate(e2(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating e2'
        endif
        allocate(e3(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating e3'
        endif
        allocate(e4(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating e4'
        endif
        allocate(e5(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating e5'
        endif
        allocate(e6(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating e6'
        endif
        allocate(e7(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating e7'
        endif
        allocate(e8(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating e8'
        endif
        allocate(f1(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating f1'
        endif
        allocate(f2(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating f2'
        endif
        allocate(f3(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating f3'
        endif
        allocate(f4(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating f4'
        endif
        allocate(f5(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating f5'
        endif
        allocate(f6(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating f6'
        endif
        allocate(f7(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating f7'
        endif
        allocate(f8(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating f8'
        endif
        allocate(bmat(9,9),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating bmat'
        endif
        allocate(bmat_coeffs(9),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating bmat_coeffs'
        endif
        allocate(bmat_evals(9),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating bmat_evals'
        endif
        allocate(ipiv(9),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating ipiv'
        endif
 

        write(6,*) '--------------------------------------------------'
        write(6,*) 'Performing RHF SCF on given starting wavefunction'
        write(6,*) 'Using DIIS for improvement of convergence'
        write(6,*) '--------------------------------------------------'

        ! construct density matrix
        dmat(:,:) = 0.0_dp
        do l1=1,norb
            do l2=1,norb
                do l3=1,(nelec/2)
                    dmat(l2,l1) = dmat(l2,l1) + (2.0_dp*coeffs(l2,l3)*&
                        &conjg(coeffs(l1,l2)))
                enddo
            enddo
        enddo

        ! fock matrix
        fmat(:,:) = 0.0_dp
        do l1=1,norb
            do l2=1,norb
                do l3=1,norb
                    do l4=1,norb
                        fmat(l2,l1) = fmat(l2,l1) + (dmat(l4,l3)*&
                            &(umat(l2,l3,l1,l4)-(0.5_dp*umat(l2,l3,l4,l1))))
                    enddo
                enddo
            enddo
        enddo

        fmat = fmat + tmat

        iter = 0
        nerror = 0
        energy = 0.0_dp
        energy_old = 0.0_dp
        std = 0.0_dp
        dmat_old = dmat

        e1(:,:) = 0.0_dp
        e2(:,:) = 0.0_dp
        e3(:,:) = 0.0_dp
        e4(:,:) = 0.0_dp
        e5(:,:) = 0.0_dp
        e6(:,:) = 0.0_dp
        e7(:,:) = 0.0_dp
        e8(:,:) = 0.0_dp
        f1(:,:) = 0.0_dp
        f2(:,:) = 0.0_dp
        f3(:,:) = 0.0_dp
        f4(:,:) = 0.0_dp
        f5(:,:) = 0.0_dp
        f6(:,:) = 0.0_dp
        f7(:,:) = 0.0_dp
        f8(:,:) = 0.0_dp


        do
            iter = iter + 1
            write(6,*) '----------------------------------------------'
            write(6,*) 'Iteration:',iter
            write(6,*) '----------------------------------------------'


   !         ! mix in part of previous density matrix to avoid oscillation
   !         ! between solution
   !         if (iter.gt.1) then!.and.(mod(iter,2).eq.0)) then

   !             dmat = (0.95_dp*dmat_old) + ((1.0_dp-0.95_dp)*dmat)
   !             ! fock matrix
   !             fmat(:,:) = 0.0_dp
   !             do l1=1,norb
   !                 do l2=1,norb
   !                     do l3=1,norb
   !                         do l4=1,norb
   !                             fmat(l2,l1) = fmat(l2,l1) + (dmat(l4,l3)*&
   !                                 &(umat(l2,l3,l1,l4)-(0.5_dp*umat(l2,l3,l4,l1))))
   !                         enddo
   !                     enddo
   !                 enddo
   !             enddo

   !             fmat = fmat + tmat

   !         endif

            ! diagonalise fmat
            write(6,*) 'Diagonalising fmat...'
            call zheev('V','U',norb,fmat,norb,evals,rscr,rlscr,rrscr,ierr)
            if (ierr.ne.0) then
                stop 'Error diagonalising fmat'
            endif

            write(6,*) 'Eigenvectors:'
            do l1=1,norb
                write(6,'(i20)') l1
                do l2=1,norb
                    write(6,*) fmat(l2,l1)
                enddo
            enddo

            write(6,*) 'Eigenvalues'
            do l1=1,norb
                write(6,*) l1,evals(l1)
            enddo

            coeffs(:,:) = fmat(:,:)
            fdiag(:) = evals(:)

            ! construct density matrix
            dmat_old = dmat
            dmat(:,:) = 0.0_dp
            do l1=1,norb
                do l2=1,norb
                    do l3=1,(nelec/2)
                        dmat(l2,l1) = dmat(l2,l1) + (2.0_dp*coeffs(l2,l3)*&
                            &conjg(coeffs(l1,l2)))
                    enddo
                enddo
            enddo

           ! if ((iter.gt.1)) then!.and.(iter.lt.5)) then
           !     dmat = (0.93_dp*dmat_old) + (0.07_dp*dmat)
           ! endif

            ! fock matrix
            fmat(:,:) = 0.0_dp
            do l1=1,norb
                do l2=1,norb
                    do l3=1,norb
                        do l4=1,norb
                            fmat(l2,l1) = fmat(l2,l1) + (dmat(l4,l3)*&
                                &(umat(l2,l3,l1,l4)-(0.5_dp*umat(l2,l3,l4,l1))))
                        enddo
                    enddo
                enddo
            enddo

            fmat = fmat + tmat

    
            ! ground state energy and difference in density matrix as 
            ! convergence criterion
            energy_old = energy
            energy = 0.0_dp
            do l1=1,norb
                do l2=1,norb
                    energy = energy + (0.5_dp*dmat(l1,l2)*&
                        &(tmat(l2,l1)+fmat(l2,l1)))
                enddo
            enddo

            ! difference in density matrices
            std = 0.0_dp
            do l1=1,norb
                do l2=1,norb
                    std = std + (abs(dmat(l2,l1)-dmat_old(l2,l1)))
                enddo
            enddo

            if (abs(aimag(energy)).gt.1e-12_dp) then
                stop 'Error: energy is complex'
            endif

            write(6,*) '-------------------------------------------------'
            write(6,*) 'The HF ground state energy of the previous and &
                & current iterations are'
            write(6,*) energy_old,energy
            write(6,*) '--------------------------------------------------'
            write(6,*) 'The summed differences in the density matrices of the &
                &previous and current iterations are'
            write(6,*) std
            write(6,*) '---------------------------------------------------'

            ! convergence
            !if (iter.eq.1) then
            if ((abs(energy-energy_old).lt.1e-4_dp).and.&
                &(abs(std).lt.1e-2_dp)) then
                write(6,*) 'Convergence achieved: stopping calculation !'
                write(6,*) '---------------------------------------------'
                exit
            endif


            ! construct error matrices needed for diis in order to improve convergence
            ! the error matrix for the ith iteration is e_i = f_i d_i s - s d_i f_i
            ! where f_i: fock matrix of ith iteration, d_i: density matrix of ith iteration

            ! in order to use successive iterations
            nerror = mod(iter,8)

            write(6,*) 'Constructing error matrices for DIIS implementation...'
            if (nerror.eq.1) then
                write(6,*) 'Constructing e1...'
                call ErrorMatrix(e1,f1,fmat,dmat)
                do l1=1,norb
                    do l2=1,norb
                        write(6,*) l1,l2,e1(l2,l1)
                    enddo
                enddo
            elseif (nerror.eq.2) then
            !elseif (nerror.eq.3) then
                write(6,*) 'Constructing e2...'
                call ErrorMatrix(e2,f2,fmat,dmat)
                do l1=1,norb
                    do l2=1,norb
                        write(6,*) l1,l2,e2(l2,l1)
                    enddo
                enddo
            elseif (nerror.eq.3) then
            !elseif (nerror.eq.5) then
                write(6,*) 'Constructing e3...'
                call ErrorMatrix(e3,f3,fmat,dmat)
                do l1=1,norb
                    do l2=1,norb
                        write(6,*) l1,l2,e3(l2,l1)
                    enddo
                enddo
            elseif (nerror.eq.4) then
            !elseif (nerror.eq.7) then
                write(6,*) 'Constructing e4...'
                call ErrorMatrix(e4,f4,fmat,dmat)
                do l1=1,norb
                    do l2=1,norb
                        write(6,*) l1,l2,e4(l2,l1)
                    enddo
                enddo
            elseif (nerror.eq.5) then
            !elseif (nerror.eq.9) then
                write(6,*) 'Constructing e5...'
                call ErrorMatrix(e5,f5,fmat,dmat)
                do l1=1,norb
                    do l2=1,norb
                        write(6,*) l1,l2,e5(l2,l1)
                    enddo
                enddo
            elseif (nerror.eq.6) then
            !elseif (nerror.eq.11) then
                write(6,*) 'Constructing e6...'
                call ErrorMatrix(e6,f6,fmat,dmat)
                do l1=1,norb
                    do l2=1,norb
                        write(6,*) l1,l2,e6(l2,l1)
                    enddo
                enddo
            elseif (nerror.eq.7) then
            !elseif (nerror.eq.13) then
                write(6,*) 'Constructing e7...'
                call ErrorMatrix(e7,f7,fmat,dmat)
                do l1=1,norb
                    do l2=1,norb
                        write(6,*) l1,l2,e7(l2,l1)
                    enddo
                enddo
            elseif (nerror.eq.0) then
            !elseif (nerror.eq.15) then
                write(6,*) 'Constructing e8...'
                call ErrorMatrix(e8,f8,fmat,dmat)
                do l1=1,norb
                    do l2=1,norb
                        write(6,*) l1,l2,e8(l2,l1)
                    enddo
                enddo
            endif

            write(6,*) '------------------------------------------------------------'

            !write(6,*) nerror
            bmat(:,:) = 0.0_dp
        !    if (iter.ge.3) then
            ! construct B matrix but only taking every second Fock matrix
             if ((nerror.eq.1).or.(nerror.eq.3).or.(nerror.eq.5).or.(nerror.eq.7).or.&
                &(nerror.eq.2).or.(nerror.eq.4).or.(nerror.eq.6).or.(nerror.eq.8)&
                &.or.(nerror.eq.0)) then
          !  if ((nerror.eq.1).or.(nerror.eq.3).or.(nerror.eq.5).or.(nerror.eq.7).or.&
          !      &(nerror.eq.9).or.(nerror.eq.11).or.(nerror.eq.13).or.(nerror.eq.15)&
          !      &.or.(nerror.eq.0)) then
                call FormBmat(bmat,e1,e1,1,1)
                call FormBmat(bmat,e1,e2,1,2)
                call FormBmat(bmat,e1,e3,1,3)
                call FormBmat(bmat,e1,e4,1,4)
                call FormBmat(bmat,e1,e5,1,5)
                call FormBmat(bmat,e1,e6,1,6)
                call FormBmat(bmat,e1,e7,1,7)
                call FormBmat(bmat,e1,e8,1,8)
                call FormBmat(bmat,e2,e2,2,2)
                call FormBmat(bmat,e2,e3,2,3)
                call FormBmat(bmat,e2,e4,2,4)
                call FormBmat(bmat,e2,e5,2,5)
                call FormBmat(bmat,e2,e6,2,6)
                call FormBmat(bmat,e2,e7,2,7)
                call FormBmat(bmat,e2,e8,2,8)
                call FormBmat(bmat,e3,e3,3,3)
                call FormBmat(bmat,e3,e4,3,4)
                call FormBmat(bmat,e3,e5,3,5)
                call FormBmat(bmat,e3,e6,3,6)
                call FormBmat(bmat,e3,e7,3,7)
                call FormBmat(bmat,e3,e8,3,8)
                call FormBmat(bmat,e4,e4,4,4)
                call FormBmat(bmat,e4,e5,4,5)
                call FormBmat(bmat,e4,e6,4,6)
                call FormBmat(bmat,e4,e7,4,7)
                call FormBmat(bmat,e4,e8,4,8)
                call FormBmat(bmat,e5,e5,5,5)
                call FormBmat(bmat,e5,e6,5,6)
                call FormBmat(bmat,e5,e7,5,7)
                call FormBmat(bmat,e5,e8,5,8)
                call FormBmat(bmat,e6,e6,6,6)
                call FormBmat(bmat,e6,e7,6,7)
                call FormBmat(bmat,e6,e8,6,8)
                call FormBmat(bmat,e7,e7,7,7)
                call FormBmat(bmat,e7,e8,7,8)
                call FormBmat(bmat,e8,e8,8,8)

                if (iter.ge.8) then
                    beforeiter = 9
                elseif (iter.lt.8) then
                    beforeiter = iter+1!(nerror+1)/2
                endif
                !write(6,*) 'itertest',beforeiter

                bmat(beforeiter,:) = -1.0_dp
                bmat(:,beforeiter) = -1.0_dp
                bmat(beforeiter,beforeiter) = 0.0_dp

                write(6,*) 'Constructing bmat...'
                do l1=1,beforeiter!9
                    do l2=1,beforeiter!9
                        write(6,*) l2,l1,bmat(l2,l1)
                    enddo
                enddo


                bmat_coeffs(:) = 0.0_dp
                bmat_coeffs(beforeiter) = -1.0_dp
                ipiv(:) = 0
                ! zgesv computes the solution to the linear equation
                ! A*X = B
                !if (iter.gt.17) then
                if (iter.gt.17) then
                    call zgesv(beforeiter,1,bmat(1:beforeiter,1:beforeiter),beforeiter,&
                        &ipiv,bmat_coeffs(1:beforeiter),beforeiter,ierr)
                    if (ierr.ne.0) then
                        write(6,*) ierr
                        stop 'Error solving system of linear equations for DIIS implementation'
                    endif

                    write(6,*) '-----------------------------------------------------------'
                    write(6,*) 'Solutions to system of linear equations of bmat...'
                    do l1=1,9
                        write(6,*) l1,bmat_coeffs(l1)
                    enddo

                    ! Construct the new fock matrix out of the previously stored ones
                    call InterpolatedFock(fmat,bmat_coeffs,f1,f2,f3,f4,f5,f6,f7,f8)

                !endif

            endif
            endif


        enddo


        ! test whether coefficients are complex or not
        ints_complex = .false.
        do l1=1,norb
            do l2=1,norb
                if (aimag(coeffs(l2,l1)).gt.1e-12_dp) then
                    ints_complex = .true.
                endif
            enddo
        enddo

        if (ints_complex) then
            write(6,*) 'Coefficients are complex'
        else
            write(6,*) 'Coefficients and hence integrals are all real'
        endif


        deallocate(fmat)
        deallocate(dmat)
        deallocate(dmat_old)
        deallocate(rscr)
        deallocate(rrscr)

        deallocate(f1)
        deallocate(e1)
        deallocate(f2)
        deallocate(e2)
        deallocate(f3)
        deallocate(e3)
        deallocate(f4)
        deallocate(e4)
        deallocate(f5)
        deallocate(e5)
        deallocate(f6)
        deallocate(e6)
        deallocate(f7)
        deallocate(e7)
        deallocate(f8)
        deallocate(e8)
        deallocate(bmat)
        deallocate(bmat_coeffs)
        deallocate(bmat_evals)
        deallocate(ipiv)

    end subroutine RHF_DIIS


    subroutine ErrorMatrix(e,f,fmat,dmat)

        ! construct error matrices needed for diis in order to improve convergence
        ! the error matrix for the ith iteration is e_i = f_i d_i s - s d_i f_i
        ! where f_i: fock matrix of ith iteration, d_i: density matrix of ith iteration

        complex(dp), allocatable, intent(inout) :: e(:,:),f(:,:)
        complex(dp), allocatable, intent(in) :: fmat(:,:),dmat(:,:)
        complex(dp), allocatable :: temp1(:,:),temp2(:,:)
        integer :: ierr,l1,l2

        allocate(temp1(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating temp1'
        endif
        allocate(temp2(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating temp2'
        endif


        ! since the transformation matrix is overwritten each iteration
        ! tmat needs to be transformed as well
        temp1(:,:) = 0.0_dp
        ! f_ij d_jk = temp1_ik
        call zgemm('N','N',norb,norb,norb,(1.0_dp,0.0_dp),fmat,norb,&
        &dmat,norb,(0.0_dp,0.0_dp),temp1,norb)
        temp2(:,:) = 0.0_dp
        ! d_ij f_jk = temp2_ik
        call zgemm('N','N',norb,norb,norb,(1.0_dp,0.0_dp),dmat,norb,&
        &fmat,norb,(0.0_dp,0.0_dp),temp2,norb)


        e = temp1 - temp2
        f = fmat

        deallocate(temp1)
        deallocate(temp2)

    end subroutine ErrorMatrix

    subroutine FormBmat(bmat,emat1,emat2,indexbmat1,indexbmat2)

        ! this subroutine is to construct the elements of the bmat
        ! b_ij = e_i . e_j = trace((B^T)A)
        ! where . is the inner product of the two matrices

        integer, intent(in) :: indexbmat1,indexbmat2
        integer :: ierr
        complex(dp), allocatable, intent(inout) :: bmat(:,:)
        complex(dp), allocatable, intent(in) :: emat1(:,:),emat2(:,:)
        complex(dp), allocatable :: temp(:,:)

        allocate(temp(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating temp'
        endif

        temp(:,:) = 0.0_dp
        call zgemm('T','N',norb,norb,norb,(1.0_dp,0.0_dp),emat1,norb,&
        &emat2,norb,(0.0_dp,0.0_dp),temp,norb)
        
        bmat(indexbmat1,indexbmat2) = 0.0_dp
        do l1=1,norb
            bmat(indexbmat1,indexbmat2) = bmat(indexbmat1,indexbmat2) &
                & + temp(l1,l1)
        enddo

        ! bmat is symmetric
        bmat(indexbmat2,indexbmat1) = bmat(indexbmat1,indexbmat2)
 
        deallocate(temp)


    end subroutine FormBmat

    subroutine InterpolatedFock(f,b,f1,f2,f3,f4,f5,f6,f7,f8)

        ! this subroutine is to construct the new fock matrix as a linear
        ! combination of the previously stored ones as suggested for DIIS

        complex(dp), allocatable, intent(inout) :: f(:,:)
        complex(dp), allocatable, intent(in) :: b(:)
        complex(dp), allocatable, intent(in) :: f1(:,:),f2(:,:),f3(:,:)
        complex(dp), allocatable, intent(in) :: f4(:,:),f5(:,:),f6(:,:),f7(:,:),f8(:,:)

        f(:,:) = 0.0_dp
        f = (b(1)*f1) + (b(2)*f2) + (b(3)*f3) + (b(4)*f4) + (b(5)*f5) +&
            & (b(6)*f6) + (b(7)*f7) + (b(8)*f8)


    end subroutine

    subroutine FindGradient(cmat,grade1,grade2)

        ! This subroutine is to find the electronic gradient in orbital-base HF theory
        ! and the electronic hessian in orbital-based HF theory
        ! E^(1)_mn = 2<CSF|[E_mn,H]|CSF> = 2(F_mn-F_nm)
        ! where F is the generalised Fock matrix
        ! F_mn = \sum_q D_mq h_nq + \sum_qrs d_mqrs g_nqrs
        ! where the singlet excitation operator is 
        ! E_pq = a^!_p\alpha a_q\alpha + a^!_p\beta a_q\beta and
        ! h_nq = <n|h|q>, g_pqrs = <pr|qs> = (pq|rs) and the 1- and 2-RDM are
        ! D_pq = <CSF|E_pq|CSF>
        ! d_pqrs = <CSF|e_pqrs|CSF> = <CSF|E_pqE_rs -\delta_rq E_ps|CSF> with
        ! e_pqrs = \sum_\sigma\tau a^!_p\sigma a^!_r\tau a_s\tau a_q\sigma (\tau.sigma
        ! are spin)

        complex(dp), allocatable, intent(in) :: cmat(:,:)
        complex(dp), allocatable :: d1mat(:,:),d2mat(:,:,:,:)
        complex(dp), allocatable :: genfmat(:,:),ymat(:,:,:,:)
        complex(dp), allocatable, intent(inout) :: grade1(:,:),grade2(:,:,:,:)
        complex(dp) :: trace_1rdm,trace_2rdm
        integer :: l1,l2,l3,l4,l5,l6,l7,l8,ierr

        allocate(d1mat(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating d1mat'
        endif
        allocate(d2mat(norb,norb,norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating d2mat'
        endif
        allocate(genfmat(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating genfmat'
        endif
        allocate(ymat(norb,norb,norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating ymat'
        endif
 

        write(6,*) '----------------------------------------------------'
        write(6,*) 'Calculating electronic gradient'
        write(6,*) '-----------------------------------------------------'

        write(6,*) 'Constructing 1- and 2-RDM...'

        ! constructe 1-RDM for spatial orbitals. the sum over spin orbitals
        ! always gives a factor of 2 for RHF
        ! 1-RDM
        d1mat(:,:) = 0.0_dp
        do l1=1,norb
            do l2=1,norb
                do l3=1,(nelec/2)
                    d1mat(l2,l1) = d1mat(l2,l1) + (2.0_dp*conjg(cmat(l2,l3))&
                        &*cmat(l1,l3))
                enddo
            enddo
        enddo

        ! 2-RDM
        d2mat(:,:,:,:) = 0.0_dp
        do l1=1,norb
            do l2=1,norb
                do l3=1,norb
                    do l4=1,norb
                        do l5=1,(nelec/2)
                            do l7=1,2  ! to account for spin
                                do l6=1,(nelec/2)
                                    do l8=1,2 ! to account for spin
                                        d2mat(l4,l3,l2,l1) = d2mat(l4,l3,l2,l1) +&
                                        &(conjg(cmat(l4,l5))*conjg(cmat(l2,l6))*&
                                        &cmat(l1,l6)*cmat(l3,l5))
                                    enddo
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo

        ! check: trace of 1-RDM should be nelec, trace 2-RDM should be 
        ! (1/2)*nelec*(nelec-1) = (1/2)*(N^2-N) = 0.5*\sum_pq d2mat(p,p,q,q) 
        ! - 0.5*\sum_p d1mat(p,p)
        ! 
        trace_1rdm = 0.0_dp
        do l1=1,norb
            trace_1rdm = trace_1rdm + d1mat(l1,l1)
        enddo
        write(6,*) 'Trace of 1-RDM:',trace_1rdm
        
         
        trace_2rdm = 0.0_dp
        do l1=1,norb
            do l2=1,norb
                trace_2rdm = trace_2rdm + (d2mat(l2,l2,l1,l1)*0.5_dp)
            enddo
            trace_2rdm = trace_2rdm - (d1mat(l1,l1)*0.5_dp)
        enddo
        write(6,*) 'Trace of 2-RDM:',trace_2rdm 


        ! the following is only true for a 1-determinantal wavefunction
        !d2mat(:,:,:,:) = 0.0_dp
        !do l1=1,norb
        !    do l2=1,norb
        !        do l3=1,norb
        !            do l4=1,norb
        !                d2mat(l4,l3,l2,l1) = d1mat(l4,l3)*d1mat(l2,l1)
        !            enddo
        !        enddo
        !    enddo
        !enddo
        !trace_2rdm = 0.0_dp
        !do l1=1,norb
        !    do l2=1,norb
        !        trace_2rdm = trace_2rdm + (d2mat(l2,l2,l1,l1)*0.5_dp)
        !    enddo
        !    trace_2rdm = trace_2rdm - (d1mat(l1,l1)*0.5_dp)
        !enddo
        !write(6,*) 'Trace of 2-RDM:',trace_2rdm 

        write(6,*) '--------------------------------------------------------'
        write(6,*) 'Constructing general fock matrix...'
      
        ! the generalised fock matrix is then
        genfmat(:,:) = 0.0_dp
        do l1=1,norb
            do l2=1,norb
                do l3=1,norb
                    genfmat(l2,l1) = genfmat(l2,l1) + (d1mat(l2,l3)*tmat(l1,l3))
                enddo
            enddo
        enddo

        do l1=1,norb
            do l2=1,norb
                do l3=1,norb
                    do l4=1,norb
                        do l5=1,norb
                            genfmat(l2,l1) = genfmat(l2,l1) + (d2mat(l2,l3,l4,l5)*&
                                &umat(l1,l4,l3,l5))
                        enddo
                    enddo
                enddo
            enddo
        enddo

        do l1=1,norb
            do l2=1,norb
                write(6,*) l1,l2,genfmat(l2,l1)
            enddo
        enddo

        write(6,*) '-------------------------------------------------------'
        write(6,*) 'Constructing gradient matrix grade1...'

        grade1(:,:) = 0.0_dp
        do l1=1,norb
            do l2=1,norb
                grade1(l2,l1) = 2.0_dp*(genfmat(l2,l1)-genfmat(l1,l2))
            enddo
        enddo

        do l1=1,norb
            do l2=1,norb
                write(6,*) l1,l2,grade1(l2,l1)
            enddo
        enddo

        ! gradient matrix should be (real ?) antisymmetric matrix
        do l1=1,norb
            do l2=l1,norb
                if (grade1(l1,l2).ne.(-grade1(l2,l1))) then
                    write(6,*) l1,l2,grade1(l1,l2),grade1(l2,l1)
                    stop 'Error: gradient matrix grade1 is not anti-symmetric'
                endif
            enddo
        enddo

        write(6,*) '--------------------------------------------------------'

        write(6,*) 'Constructing y matrix...'

        ! y_pqrs = \sum_mn ((d_pmrn+d_pmnr)*g_qmns + (d_prmn*g_qsmn))
        ymat(:,:,:,:) = 0.0_dp
        do l1=1,norb
            do l2=1,norb
                do l3=1,norb
                    do l4=1,norb
                        do l5=1,norb
                            do l6=1,norb
                                ymat(l4,l3,l2,l1) = ymat(l4,l3,l2,l1) + &
                                    &((d2mat(l4,l5,l2,l6) + d2mat(l4,l5,l6,l2))&
                                    &*umat(l3,l6,l5,l1)) + &
                                    &(d2mat(l4,l2,l5,l6)*umat(l3,l5,l1,l6))
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo

        write(6,*) 'Constructing hessian matrix grade2...'
        ! e^2_pqrs = (1-P_pq)*(1-P_rs)*[(2*D_prh_qs - (F_pr + F_rp)\delta_qs +
        ! 2*Y_pqrs]
        grade2(:,:,:,:) = 0.0_dp
        do l1=1,norb
            do l2=1,norb
                do l3=1,norb
                    do l4=1,norb
                        ! no permutation
                        grade2(l4,l3,l2,l1) = grade2(l4,l3,l2,l1) + (2.0_dp*d1mat(l4,l2)&
                            &*tmat(l3,l1)) + (2.0_dp*ymat(l4,l3,l2,l1))
                        ! permute l4 and l3 
                        grade2(l4,l3,l2,l1) = grade2(l4,l3,l2,l1) + (2.0_dp*d1mat(l3,l2)&
                            &*tmat(l4,l1)) + (2.0_dp*ymat(l3,l4,l2,l1))
                        ! permute l2 and l1 
                        grade2(l4,l3,l2,l1) = grade2(l4,l3,l2,l1) + (2.0_dp*d1mat(l4,l1)&
                            &*tmat(l3,l2)) + (2.0_dp*ymat(l4,l3,l1,l2))
                        ! permute l4 and l3,as well as l2 and l1 
                        grade2(l4,l3,l2,l1) = grade2(l4,l3,l2,l1) + (2.0_dp*d1mat(l3,l1)&
                            &*tmat(l4,l2)) + (2.0_dp*ymat(l3,l4,l1,l2))
                        if (l3.eq.l1) then
                            ! no permutation
                            grade2(l4,l3,l2,l1) = grade2(l4,l3,l2,l1) &
                                &- (genfmat(l4,l2)+genfmat(l2,l4))
                        endif
                        if (l4.eq.l1) then
                            ! permute l4 and l3
                            grade2(l4,l3,l2,l1) = grade2(l4,l3,l2,l1) &
                                &- (genfmat(l3,l2)+genfmat(l2,l3))
                        endif
                        if (l3.eq.l2) then
                            ! permute l2 and l1
                            grade2(l4,l3,l2,l1) = grade2(l4,l3,l2,l1) &
                                &- (genfmat(l4,l1)+genfmat(l1,l4))
                        endif
                        if (l4.eq.l2) then
                            ! permute l4 and l3, as well as l2 and l1
                            grade2(l4,l3,l2,l1) = grade2(l4,l3,l2,l1) &
                                &- (genfmat(l3,l1)+genfmat(l1,l3))
                        endif
                    enddo
                enddo
            enddo
        enddo

        do l1=1,norb
            do l2=1,norb
                do l3=1,norb
                    do l4=1,norb
                        write(6,*) l1,l2,l3,l4,grade2(l4,l3,l2,l1)
                    enddo
                enddo
            enddo
        enddo

        write(6,*) '--------------------------------------------------------'
 

        deallocate(d1mat)
        deallocate(d2mat)
        deallocate(genfmat)
        !deallocate(grade1)
        deallocate(ymat)
        !deallocate(grade2)

    end subroutine FindGradient

    subroutine TransformFromGradient(cmat,g1,alpha)

        ! this subroutine is to generate a special orthogonal matrix R which performs
        ! a unitary transformation of the orbitals using the gradient, such that
        ! the transformation follows the steepest descent path in order to find 
        ! another HF energy minimum
        ! transformation R = exp(rX) where rX is a real-antisymmetric matrix
        ! diagonalisation of rX gives rX=iV \tau V^! where \tau is a diagonal matrix

        complex(dp), allocatable, intent(in) :: g1(:,:)
        complex(dp), allocatable, intent(inout) :: cmat(:,:)
        complex(dp), allocatable :: vmat(:,:)
        real(dp), allocatable :: squarexmat(:,:),wmat(:,:),sxevals(:)
        real(dp), allocatable :: sintau(:,:),costau(:,:)
        real(dp), allocatable :: rrscr(:)
        real(dp), allocatable :: rmat(:,:),rmat1(:,:),rmat2(:,:)
        real(dp), allocatable :: temp1(:,:),temp2(:,:)
        real(dp), intent(in) :: alpha
        integer :: l1,l2,l3,l4,ierr,rlscr


        ! when X is anti-hermitian, iX is hermitian (needed in order to use zheev)
        ! from which it follows iX = V \tau V^! and X = V (-i\tau) V^!

        rlscr = 4*norb
        allocate(vmat(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating vmat'
        endif
        allocate(rrscr(rlscr),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating rrscr'
        endif
        allocate(squarexmat(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating squarexmat'
        endif
        allocate(wmat(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating wmat'
        endif
        allocate(sxevals(norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating sxevals'
        endif
        allocate(sintau(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating sintau'
        endif
        allocate(costau(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating costau'
        endif
        allocate(rmat(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating rmat'
        endif
        allocate(rmat1(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating rmat1'
        endif
        allocate(rmat2(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating rmat2'
        endif
        allocate(temp1(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating temp1'
        endif
        allocate(temp2(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating temp2'
        endif
      
        write(6,*) '------------------------------------------------------------'
        write(6,*) 'Constructing transformation matrix U=exp(-X) for orbital rotation'
        write(6,*) 'using electronic HF energy gradient as X assuming skew-symmetric X'
        write(6,*) '------------------------------------------------------------'

        ! in order to avoid complex arithmetic
        ! X^2 is used which is symmetric and therefore has 
        ! real eigenvalues and eigenvectors

        ! need to form RX from X since X only contains real numbers
        ! however, X is anti-Hermitian and RX needs to be skew-symmetric (anti-symmetric)
        ! hence RX = (X - X^T)/2 is needed

        vmat(:,:) = 0.0_dp
        vmat = g1 - transpose(g1)
        vmat = vmat/2.0_dp

        wmat = -alpha*real(vmat,dp)
        squarexmat(:,:) = 0.0_dp
        ! X^2: wmat_ij wmat_jk = wmat_ik
        call dgemm('N','N',norb,norb,norb,1.0_dp,wmat,norb,&
        &wmat,norb,0.0_dp,squarexmat,norb)

        ! symmetric matrix: A = A^T
        write(6,*) 'Matrix X^2:'
        do l1=1,norb
            do l2=1,norb
                write(6,'(2i5,1G20.12)') l2,l1,squarexmat(l2,l1)
            enddo
        enddo

        write(6,*) '------------------------------------------------------------------'

        call dsyev('V','U',norb,squarexmat,norb,sxevals,rrscr,rlscr,ierr)
        if (ierr.ne.0) then
            stop 'Error diagonalising X^2'
        endif

        write(6,*) 'Eigenvectors of X^2:'
        do l1=1,norb
            write(6,'(i9)') l1
            do l2=1,norb
                write(6,'(i5,G20.12)') l2,squarexmat(l2,l1)
            enddo
        enddo
        write(6,*) 'Eigenvalues of X^2:'
        do l1=1,norb
            write(6,'(i5,G20.12)') l1,sxevals(l1)
        enddo


        write(6,*) '--------------------------------------------------------'
        ! cos(tau)
        costau(:,:) = 0.0_dp
        do l1=1,norb
            costau(l1,l1) = cos(sqrt(abs(sxevals(l1))))
        enddo

        ! sin(tau) * (tau^-1)
        sintau(:,:) = 0.0_dp
        do l1=1,norb
            if (abs(sxevals(l1)).gt.1e-8_dp) then
               sintau(l1,l1) = sin(sqrt(abs(sxevals(l1))))/sqrt(abs(sxevals(l1)))
            else
                ! if tau = 0: sin(tau)*(tau^-1) = 1
                sintau(l1,l1) = 1.0_dp
            endif
        enddo


        write(6,*) 'Diagonal elements of cos(tau) and sin(tau) (tau_i is ith eigenvalue of X^2)'
        do l1=1,norb
            write(6,'(i5,2G20.12)') l1,costau(l1,l1),sintau(l1,l1)
        enddo

        write(6,*) '------------------------------------------------------------------'
        write(6,*) 'Transformation matrix R=exp(-X)'

        ! R = W cos(tau) W^T + W tau^-1 sin(tau) W^T X
        ! where X^2 = -W tau^2 W^T

        temp1(:,:) = 0.0_dp
        ! cos(tau)_ij (W_kj)^T = temp1_ik
        call dgemm('N','T',norb,norb,norb,1.0_dp,costau,norb,&
        &squarexmat,norb,0.0_dp,temp1,norb)
        ! W_ij temp1_jk = rmat_ik
        rmat1(:,:) = 0.0_dp
        call dgemm('N','N',norb,norb,norb,1.0_dp,squarexmat,norb,&
        &temp1,norb,0.0_dp,rmat1,norb)


        temp1(:,:) = 0.0_dp
        ! (W_ij)^T X_ik = temp1_jk
        temp1(:,:) = 0.0_dp
        call dgemm('T','N',norb,norb,norb,1.0_dp,squarexmat,norb,&
        &wmat,norb,0.0_dp,temp1,norb)
        ! sintau_ij temp1_jk = temp2_ik
        temp2(:,:) = 0.0_dp
        call dgemm('N','N',norb,norb,norb,1.0_dp,sintau,norb,&
        &temp1,norb,0.0_dp,temp2,norb)
        ! X_ij temp1_jk = rmat2_ik
        rmat2(:,:) = 0.0_dp
        call dgemm('N','N',norb,norb,norb,1.0_dp,squarexmat,norb,&
        &temp2,norb,0.0_dp,rmat2,norb)


        rmat = rmat1 + rmat2
        
        do l1=1,norb
            do l2=1,norb
                write(6,'(2i5,G20.12)') l2,l1,rmat(l2,l1)
            enddo
        enddo

        ! test for unitarity
        call dgemm('N','T',norb,norb,norb,1.0_dp,rmat,norb,&
        &rmat,norb,0.0_dp,temp2,norb)

        write(6,*) '---------------------------------------------------------'
        write(6,*) 'Test for R being unitary'
        do l1=1,norb
            do l2=1,norb
                if (abs(temp2(l2,l1)).gt.1e-12_dp) then
                    write(6,*) l2,l1,temp2(l2,l1)
                endif
            enddo
        enddo

        write(6,*) '----------------------------------------------------------'

        ! return transformation matrix
        cmat(:,:) = 0.0_dp
        cmat = rmat

        deallocate(rrscr)
        deallocate(squarexmat)
        deallocate(wmat)
        deallocate(vmat)
        deallocate(sxevals)
        deallocate(sintau)
        deallocate(costau)
        deallocate(temp1)
        deallocate(temp2)
        deallocate(rmat1)
        deallocate(rmat2)
        deallocate(rmat)

    end subroutine TransformFromGradient

    subroutine TransformFromGradient_Complex(cmat,g1,alpha)

        ! This subroutine is simialr to TransformFromGradient with then only
        ! difference that it uses complex arithmetics involving a more general
        ! approach
        ! This subroutine is to generate a special orthogonal matrix R which performs
        ! a unitary transformation of the orbitals using the gradient, such that
        ! the transformation follows the steepest descent path in order to find 
        ! another HF energy minimum
        ! transformation R = exp(rX) where rX is a real-antisymmetric matrix
        ! diagonalisation of rX gives rX=iV \tau V^! where \tau is a diagonal matrix

        complex(dp), allocatable, intent(in) :: g1(:,:)
        complex(dp), allocatable, intent(inout) :: cmat(:,:)
        complex(dp), allocatable :: wmat(:,:),diagexp(:,:),temp1(:,:)
        real(dp), allocatable :: wevals(:)
        real(dp), allocatable :: rrscr(:),rrlscr(:)
        complex(dp), allocatable :: rmat(:,:),rmat1(:,:)
        real(dp), intent(in) :: alpha
        integer :: l1,l2,l3,l4,ierr,rlscr


        ! when X is anti-hermitian, iX is hermitian (needed in order to use zheev)
        ! from which it follows iX = V \tau V^! and X = V (-i\tau) V^!

        rlscr = 4*norb
        allocate(wmat(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating wmat'
        endif
        allocate(rrscr(rlscr),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating rrscr'
        endif
        allocate(diagexp(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating diagexp'
        endif
        allocate(wevals(norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating wevals'
        endif
        allocate(rmat(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating rmat'
        endif
        allocate(rmat1(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating rmat1'
        endif
        allocate(temp1(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating temp1'
        endif
        allocate(rrlscr(rlscr),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating rlrscr'
        endif
      
        write(6,*) '------------------------------------------------------------'
        write(6,*) 'Constructing transformation matrix U=exp(-X) for orbital rotation'
        write(6,*) 'using electronic HF energy gradient as X assumin anti-Hermitian X'
        write(6,*) '------------------------------------------------------------'

        ! the matrix exponential R = exp(x)
        ! a unitary matrix X can be diagonalised using
        ! X = V w V^!
        ! X is anti-Hermitian X_ij = - (X_ji)^!
        ! where w^! is the hermitian conjugate and 
        ! w is a diagonal matrix with i delta_k where i delta_k are 
        ! the purely imaginary eigenvalues of the anti-Hermitian matrix X
        ! the matrix expontian is then 
        ! exp(X) = V exp(i delta) V^!
        ! in order to diagonalise X more easily iX is diagonalised
        ! which is Hermitian and the matrix expontial is then
        ! X = V exp(-i delta) V^!

        ! Hermitian iX
        wmat = -alpha*imaginaryi*g1

        ! Hermitain iX
        write(6,*) 'Hermitian matrix iX:'
        do l1=1,norb
            do l2=1,norb
                write(6,'(2i5,2G20.12)') l2,l1,wmat(l2,l1)
            enddo
        enddo

        write(6,*) '------------------------------------------------------------------'

        call zheev('V','U',norb,wmat,norb,wevals,rrscr,rlscr,rrlscr,ierr)
        if (ierr.ne.0) then
            stop 'Error diagonalising iX'
        endif
 
        write(6,*) 'Eigenvectors of iX:'
        do l1=1,norb
            write(6,'(i9)') l1
            do l2=1,norb
                write(6,'(i4,2G20.12)') l2,wmat(l2,l1)
            enddo
        enddo
        write(6,*) 'Eigenvalues of iX:'
        do l1=1,norb
            write(6,'(i5,G20.12)') l1,wevals(l1)
        enddo

        ! exp(-i delta)
        diagexp(:,:) = 0.0_dp
        do l1=1,norb
            diagexp(l1,l1) = exp(-imaginaryi*wevals(l1))
        enddo

        write(6,*) '---------------------------------------------------------'
        write(6,*) 'Diagonal elements of exp(i delta)'
        do l1=1,norb
            write(6,'(i5,2G20.12)') l1,diagexp(l1,l1)
        enddo

        write(6,*) '------------------------------------------------------------------'
        write(6,*) 'Transformation matrix R=exp(-X)'
       
        ! V exp(-idelta) V^!
        rmat1(:,:) = 0.0_dp
        call zgemm('N','C',norb,norb,norb,(1.0_dp,0.0_dp),diagexp,norb,&
        &wmat,norb,(0.0_dp,0.0_dp),rmat1,norb)
        rmat(:,:) = 0.0_dp
        call zgemm('N','N',norb,norb,norb,(1.0_dp,0.0_dp),wmat,norb,&
        &rmat1,norb,(0.0_dp,0.0_dp),rmat,norb)


        do l1=1,norb
            do l2=1,norb
                write(6,'(2i5,2G20.12)') l2,l1,rmat(l2,l1)
            enddo
        enddo

        ! test for imaginary contributions
        do l1=1,norb
            do l2=1,norb
                if (abs(aimag(rmat(l2,l1))).gt.1e-12_dp) then
                    write(6,*) 'WARNING: rotation matrix contains complex elements'
                endif
            enddo
        enddo

        temp1(:,:) = 0.0_dp
        call zgemm('N','C',norb,norb,norb,(1.0_dp,0.0_dp),rmat,norb,&
        &rmat,norb,(0.0_dp,0.0_dp),temp1,norb)


        write(6,*) '---------------------------------------------------------'
        write(6,*) 'Test for R being unitary'
        do l1=1,norb
            do l2=1,norb
                if (abs(temp1(l2,l1)).gt.1e-12_dp) then
                    write(6,*) l2,l1,temp1(l2,l1)
                endif
            enddo
        enddo

        write(6,*) '----------------------------------------------------------'


        ! return transformation matrix
        cmat(:,:) = 0.0_dp
        cmat = rmat

        deallocate(rrscr)
        deallocate(rrlscr)
        deallocate(diagexp)
        deallocate(wmat)
        deallocate(wevals)
        deallocate(rmat1)
        deallocate(temp1)
        deallocate(rmat)

    end subroutine TransformFromGradient_Complex


    subroutine SteepestDescent(cmat,grade1,grade2)

        ! this subroutine is to perform a similar operation as steepest destcent
        ! in order to go downhill to a new minimum on the HF energy surface if
        ! an instability in the HF energy has been found

        complex(dp), allocatable, intent(inout) :: cmat(:,:)
        complex(dp), allocatable, intent(inout) :: grade1(:,:),grade2(:,:,:,:)
        complex(dp), allocatable :: amat(:,:),fmat(:,:),dmat(:,:),energies(:)
        !complex(dp), allocatable :: hfdmat(:,:)
        complex(dp), allocatable :: hfcoeffs(:,:)
        real(dp), allocatable :: nos(:,:),no_occ(:),fevals(:)
        real(dp), allocatable :: norrscr(:),rrscr(:),rrlscr(:)
        complex(dp) :: energy,energy_old!,normalization
        real(dp) :: cauchyschwarz!,diffcauchyschwarz
        real(dp) :: minenergy
        integer :: minenergyloc(1)
        real(dp), allocatable :: alpha(:)
        integer :: l1,l2,l3,l4,l5,ierr,rlscr

        rlscr = 4*norb
        allocate(amat(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating amat'
        endif
       ! allocate(fmat(norb,norb),stat=ierr)
       ! if (ierr.ne.0) then
       !     stop 'Error allocating fmat'
       ! endif
        allocate(dmat(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating dmat'
        endif
        allocate(alpha(200),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating alpha'
        endif
        allocate(energies(200),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating energies'
        endif
        allocate(hfcoeffs(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating hfcoeffs'
        endif
        allocate(norrscr(rlscr),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating norrscr'
        endif
        allocate(nos(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating nos'
        endif
        allocate(no_occ(norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating no_occ'
        endif
     !   allocate(fevals(norb),stat=ierr)
     !   if (ierr.ne.0) then
     !       stop 'Error allocating fevals'
     !   endif
 

        write(6,*) '-------------------------------------------------------------------'
        write(6,*) 'Approaching new HF energy minimum by following gradient'
        write(6,*) '-------------------------------------------------------------------'


        ! the vector corresponding to the lowest eigenvalue of the instability
        ! matrix give the direction of steepest descent, i.e. the gradient in
        ! this direction

        write(6,*) 'Performing orbital rotation following the direction of the gradient'

        ! the original HF orbitals as read in from the FCIDUMP file serve as
        ! basis set
        hfcoeffs(:,:) = 0.0_dp
        do l1=1,norb!(nelec/2)
            hfcoeffs(l1,l1) = 1.0_dp
        enddo

        grade1(:,:) = 0.0_dp
        grade1 = cmat

     !   write(6,*) 'Xmatrix for tranformation'
     !   do l1=1,norb
     !       do l2=1,norb
     !           write(6,*) l2,l1,grade1(l2,l1)
     !       enddo
     !   enddo
        ! values of alpha to consider in R = exp(-alpha X)
        alpha(:) = 0.0_dp
        alpha(1) = -5.0_dp
        ! alpha values around 0 are not considered
        do l1=2,100
            alpha(l1) = alpha(l1-1) + 0.04_dp
        enddo
        alpha(101) = 1.0_dp
        do l1=102,200
            alpha(l1) = alpha(l1-1) + 0.04_dp
        enddo


        ! calculate energy as a function of alpha
        do l5=1,200

            write(6,*) 'Value for parameter alpha:',alpha(l5)
        
 
            !call TransformFromGradient(cmat,grade1,1.0_dp)
            call TransformFromGradient_Complex(cmat,grade1,alpha(l5))

            ! the new orbitals are found by transforming the old orbitals
            ! this amount to a matrix multplication in order to obtain the 
            ! new coefficients 
            ! phi_i = \sum_j R_ij \sum_\mu C_j\mu psi_\mu = \sum_j\sum R_ij C_j\mu psi_\mu
            ! and the new coefficients are A_ik = R_ij C_jk

            ! the ith column in cmat corresponds to the transformation of the ith hf orbital
            amat(:,:) = 0.0_dp
            ! a_ij = (r_ji)^T C_jk
            call zgemm('T','N',norb,norb,norb,(1.0_dp,0.0_dp),cmat,norb,&
            &hfcoeffs,norb,(0.0_dp,0.0_dp),amat,norb)

            amat = cmat

            ! tansformed orbitals
            write(6,*) 'The transformed orbitals are'
            do l1=1,norb
                write(6,'(i9)') l1
                do l2=1,norb
                    write(6,'(i5,2g20.12)') l2,amat(l2,l1)
                enddo
            enddo

            write(6,*) '---------------------------------------------------------'

            ! construct density and fock matrices in order to evaluate the energy

            ! construct density matrix
            dmat(:,:) = 0.0_dp
            do l1=1,norb
                do l2=1,norb
                    do l3=1,(nelec/2)
                        dmat(l2,l1) = dmat(l2,l1) + (2.0_dp*amat(l1,l3)*&
                            &conjg(amat(l2,l3)))
                    enddo
                enddo
            enddo

            write(6,*) 'Density matrix for transformed orbitals'
            do l1=1,norb
                do l2=1,norb
                    write(6,*) l2,l1,dmat(l2,l1)
                enddo
            enddo

            ! if the fock matrix should be evaluated at diagonalised at this point

   !     ! fock matrix
   !     fmat(:,:) = 0.0_dp
   !     do l1=1,norb
   !         do l2=1,norb
   !             do l3=1,norb
   !                 do l4=1,norb
   !                     fmat(l2,l1) = fmat(l2,l1) + (dmat(l3,l4)*&
   !                         &(umat(l2,l3,l1,l4)-(0.5_dp*umat(l2,l3,l4,l1))))
   !                 enddo
   !             enddo
   !         enddo
   !     enddo

   !     fmat = fmat + tmat
   !     call zheev('V','U',norb,fmat,norb,fevals,rrscr,rlscr,rrlscr,ierr)
   !     if (ierr.ne.0) then
   !         stop 'Error diagonalising fmat'
   !     endif

   !     write(6,*) 'Eigenvectors of fock matrix:'
   !     do l1=1,norb
   !         write(6,*) l1
   !         do l2=1,norb
   !             write(6,*) l2,fmat(l2,l1)
   !         enddo
   !     enddo

   !     write(6,*) 'Eigenvalues of fock matrix'
   !     do l1=1,norb
   !         write(6,*) l1,fevals(l1)
   !     enddo
 
   !     ! construct density matrix
   !     hfdmat(:,:) = 0.0_dp
   !     do l1=1,norb
   !         do l2=1,norb
   !             do l3=1,(nelec/2)
   !                 hfdmat(l2,l1) = hfdmat(l2,l1) + (2.0_dp*fmat(l1,l3)*&
   !                     &conjg(fmat(l2,l3)))
   !             enddo
   !         enddo
   !     enddo


            ! test that the new density matrix fulfills the Cauchy-Schwarz inequality
            ! |d_pq| =< (|d_pp||d_qq|)^(-1/2)
           do l1=1,norb
                do l2=1,norb
                    cauchyschwarz = sqrt(abs(real(dmat(l1,l1),dp))*abs(real(dmat(l2,l2),dp)))
                    if (abs(cauchyschwarz).lt.abs(real(dmat(l2,l1),dp))) then
                        write(6,*) 'Cauchy-Schwarz is violated:',dmat(l2,l1),dmat(l2,l2),&
                        &dmat(l1,l1),cauchyschwarz
                    endif
                enddo
            enddo

            ! find natural orbitals
            nos = real(dmat,dp)
            call dsyev('V','U',norb,nos,norb,no_occ,norrscr,rlscr,ierr)
            if (ierr.ne.0) then
                stop 'Error diagonalising nos'
            endif
     
            write(6,*) 'Natural Orbitals:'
            do l1=1,norb
                write(6,*) l1
                do l2=1,norb
                    write(6,*) l2,nos(l2,l1)
                enddo
            enddo

            write(6,*) 'Occupation numbers:'
            do l1=1,norb
                write(6,*) l1,no_occ(l1)
            enddo

            write(6,*) 'Sum of natural orbital occupation numbers:',sum(no_occ)

            ! calculate energy
            energy = 0.0_dp
            do l1=1,norb
                do l2=1,norb
                    do l3=1,norb
                        do l4=1,norb
                            ! Coulomb energy
                            energy = energy + (dmat(l4,l2)*dmat(l3,l1)*umat(l4,l3,l2,l1))
                            ! Exchange energy
                            energy = energy - (0.5_dp*dmat(l4,l2)*&
                                    &dmat(l3,l1)*umat(l4,l3,l1,l2))
                        enddo
                    enddo
                enddo
            enddo
            energy = 0.5_dp*energy

            ! kinetic energy
            do l1=1,norb
                do l2=1,norb
                    energy = energy + (dmat(l2,l1)*tmat(l2,l1))
                enddo
            enddo


            write(6,*) 'The new HF energy is:',energy

            energies(l5) = energy

        enddo
        
        write(6,*) '---------------------------------------------------------'
        write(6,*) 'The HF energy E(alpha) is: (alpha,E(alpha))'
        do l1=1,200
            write(6,*) l1,alpha(l1),energies(l1)
        enddo

        ! choose the alpha which gives the lowest energy
        minenergy = minval(real(energies,dp))
        minenergyloc = minloc(real(energies,dp))
        

        write(6,*) '-------------------------------------------------------------'
        write(6,*) 'Choosing the transformation giving the energy:',minenergyloc,minenergy
        ! comput transformation
        call TransformFromGradient(cmat,grade1,alpha(minenergyloc(1)))

        amat = cmat
        ! values of alpha to consider in R = exp(-alpha X)
        alpha(:) = 0.0_dp
        alpha(1) = -2.0_dp
        do l1=2,200
            alpha(l1) = alpha(l1-1) + 0.02_dp
        enddo

        energy_old = energy

        ! find gradient at this new point
        write(6,*) 'Evaluating electronic gradient at this point...'
        call FindGradient(amat,grade1,grade2)

        ! update coefficients of basis set
        hfcoeffs = amat

        ! proceed via steepest descent until the change in the gradient is small
        ! and either a diis or a newton-raphson approach is better

        do

            ! calculate energy as a function of alpha
            do l5=1,200

                write(6,*) 'Value for parameter alpha',alpha(l5)
            
                call TransformFromGradient_Complex(cmat,grade1,alpha(l5))
               
               ! since the transform one electron integrals
                ! tmat needs to be transformed as well
                amat(:,:) = 0.0_dp
                ! a_ij = r_ij C)jk
                call zgemm('T','N',norb,norb,norb,(1.0_dp,0.0_dp),cmat,norb,&
                &hfcoeffs,norb,(0.0_dp,0.0_dp),amat,norb)
                

                ! construct density matrix
                dmat(:,:) = 0.0_dp
                do l1=1,norb
                    do l2=1,norb
                        do l3=1,(nelec/2)
                            dmat(l2,l1) = dmat(l2,l1) + (2.0_dp*amat(l1,l3)*&
                                &conjg(amat(l2,l3)))
                        enddo
                    enddo
                enddo

                write(6,*) 'Density matrix for transformed orbitals'
                do l1=1,norb
                    do l2=1,norb
                        write(6,*) l2,l1,dmat(l2,l1)
                    enddo
                enddo

                ! if the fock matrix is to be evaluated at this point

       !         ! fock matrix
       !         fmat(:,:) = 0.0_dp
       !         do l1=1,norb
       !             do l2=1,norb
       !                 do l3=1,norb
   !                     do l4=1,norb
   !                         fmat(l2,l1) = fmat(l2,l1) + (dmat(l3,l4)*&
   !                             &(umat(l2,l3,l1,l4)-(0.5_dp*umat(l2,l3,l4,l1))))
   !                     enddo
   !                 enddo
   !             enddo
   !         enddo

   !         fmat = fmat + tmat
   !         call zheev('V','U',norb,fmat,norb,fevals,rrscr,rlscr,rrlscr,ierr)
   !         if (ierr.ne.0) then
   !             stop 'Error diagonalising fmat'
   !         endif

   !         write(6,*) 'Eigenvectors of fock matrix:'
   !         do l1=1,norb
   !             write(6,*) l1
   !             do l2=1,norb
   !                 write(6,*) l2,fmat(l2,l1)
   !             enddo
   !         enddo

   !         write(6,*) 'Eigenvalues of fock matrix'
   !         do l1=1,norb
   !             write(6,*) l1,fevals(l1)
   !         enddo
   !  
   !         ! construct density matrix
   !         hfdmat(:,:) = 0.0_dp
   !         do l1=1,norb
   !             do l2=1,norb
   !                 do l3=1,(nelec/2)
   !                     hfdmat(l2,l1) = hfdmat(l2,l1) + (2.0_dp*fmat(l1,l3)*&
   !                         &conjg(fmat(l2,l3)))
   !                 enddo
   !             enddo
   !         enddo


                ! test that density matrix fulfills Cauchy-Schwarz inequality
                ! |d_pq| =< (|d_pp||d_qq|)^(-1/2)
     
                do l1=1,norb
                    do l2=1,norb
                        cauchyschwarz = sqrt(abs(real(dmat(l1,l1),dp))*abs(real(dmat(l2,l2),dp)))
                        if (abs(cauchyschwarz).lt.abs(real(dmat(l2,l1),dp))) then
                            write(6,*) 'Cauchy-Schwarz is violated!',dmat(l2,l1),dmat(l2,l2),&
                            &dmat(l1,l1),cauchyschwarz
                        endif
                    enddo
                enddo

                ! find natural orbitals
                nos = real(dmat,dp)
                call dsyev('V','U',norb,nos,norb,no_occ,norrscr,rlscr,ierr)
                if (ierr.ne.0) then
                    stop 'Error diagonalising nos'
                endif

                write(6,*) 'Natural Orbitals:'
                do l1=1,norb
                    write(6,*) l1
                    do l2=1,norb
                        write(6,*) l2,nos(l2,l1)
                    enddo
                enddo

                write(6,*) 'Occupation numbers:'
                do l1=1,norb
                    write(6,*) l1,no_occ(l1)
                enddo

                write(6,*) 'Sum of natural orbital occupation numbers:',sum(no_occ)

                ! calculate energy   
                energies(l5) = 0.0_dp
                do l1=1,norb
                    do l2=1,norb
                        do l3=1,norb
                            do l4=1,norb
                                ! Coulomb energy
                                energies(l5) = energies(l5) + (dmat(l4,l2)*dmat(l3,l1)*umat(l4,l3,l2,l1))
                                ! exchange energy
                                energies(l5) = energies(l5) - (0.5_dp*dmat(l4,l2)*&
                                    &dmat(l3,l1)*umat(l4,l3,l1,l2))
                            enddo
                        enddo
                    enddo
                enddo
                energies(l5) = 0.5_dp*energies(l5)

                ! kinetic energy
                do l1=1,norb
                    do l2=1,norb
                        energies(l5) = energies(l5) + (dmat(l2,l1)*tmat(l2,l1))
                    enddo
                enddo

            enddo

            write(6,*) '--------------------------------------------------------------'

            write(6,*) 'The HF energy E(alpha) is: (alpha,E(alpha))'
            do l1=1,200
                write(6,*) l1,alpha(l1),energies(l1)
            enddo

            ! choose the alpha which gives the lowest energy
            minenergy = minval(real(energies,dp))
            minenergyloc = minloc(real(energies,dp))
            

            ! only go downhill
            if (minenergy.gt.real(energy_old,dp)) then
                write(6,*) 'Energy increases: stopping steepest-descent method'
                exit
            endif

            write(6,*) '-------------------------------------------------------------'
            write(6,*) 'Choosing the transformation giving the energy:',minenergyloc,minenergy
            ! comput transformation
            call TransformFromGradient(cmat,grade1,alpha(minenergyloc(1)))
        
            ! since the transform one electron integrals
            ! tmat needs to be transformed as well
            amat(:,:) = 0.0_dp
            ! a_ij = r_ij C)jk
            call zgemm('T','N',norb,norb,norb,(1.0_dp,0.0_dp),cmat,norb,&
            &hfcoeffs,norb,(0.0_dp,0.0_dp),amat,norb)

            hfcoeffs = amat
             ! construct density matrix
            dmat(:,:) = 0.0_dp
            do l1=1,norb
                do l2=1,norb
                    do l3=1,(nelec/2)
                        dmat(l2,l1) = dmat(l2,l1) + (2.0_dp*amat(l1,l3)*&
                            &conjg(amat(l2,l3)))
                    enddo
                enddo
            enddo


            write(6,*) 'Density matrix for transformed orbitals'
            do l1=1,norb
                do l2=1,norb
                    write(6,*) l2,l1,dmat(l2,l1)
                enddo
            enddo
            ! find natural orbitals
            nos = real(dmat,dp)
            call dsyev('V','U',norb,nos,norb,no_occ,norrscr,rlscr,ierr)
            if (ierr.ne.0) then
                stop 'Error diagonalising nos'
            endif

            write(6,*) 'Natural Orbitals:'
            do l1=1,norb
                write(6,*) l1
                do l2=1,norb
                    write(6,*) l2,nos(l2,l1)
                enddo
            enddo

            write(6,*) 'Occupation numbers:'
            do l1=1,norb
                write(6,*) l1,no_occ(l1)
            enddo

            write(6,*) 'Sum of natural orbital occupation numbers:',sum(no_occ)


            ! decrease in energy is only small and continuation with another method would be
            ! better
            if (abs(energy_old - minenergy).lt.1e-5_dp) then
                write(6,*) 'Decrease in energy is only small: stopping &
                    &steepest-descent method'
                exit
            endif

            energy_old = minenergy

        enddo

        coeffs = real(amat,dp)

        ! final coefficients
        write(6,*) '--------------------------------------------------------------------'
        write(6,*) 'The final coefficients:'

        do l1=1,norb
            do l2=1,norb
                write(6,*) l2,l1,coeffs(l2,l1)
            enddo
        enddo

        deallocate(amat)
        !deallocate(fmat)
        !deallocate(fevals)
        deallocate(nos)
        deallocate(hfcoeffs)
        deallocate(no_occ)
        deallocate(dmat)
        deallocate(alpha)
        deallocate(energies)
        deallocate(norrscr)

    end subroutine SteepestDescent


    subroutine NewtonRaphson(cmat,g1,g2)

        ! This subroutine is to apply a Newton-Raphson approach in order to find a 
        ! minimum
        ! using X = -(H-fI)^-1 g
        ! where g is the gradient, H the hessian, I the identity and f and shift
        ! parameter in order to deal with zero eigenvalues of the Hessian

        complex(dp), allocatable, intent(in) :: g1(:,:),g2(:,:,:,:)
        complex(dp), allocatable, intent(inout) :: cmat(:,:)
        real(dp), allocatable :: hessian(:,:),gradient(:),xmat(:)
        real(dp), allocatable :: rrscr(:),hessian_evals(:),hessian_eigenvec(:,:)
        real(dp), allocatable :: hessian_coeffs(:),xfromh(:)
        real(dp), allocatable :: inversehessian(:,:),test(:,:)
        integer :: l1,l2,l3,l4,ierr,m,n,nexcit,lrscr
        real(dp) :: mineval

        ! number of matrix elements
        nexcit = norb*norb
        lrscr = 4*nexcit

        write(6,*) '-----------------------------------------------------------'
        write(6,*) 'Newton-Raphson Step for finding a minima in the HF energy'
        write(6,*) '-----------------------------------------------------------'

        allocate(hessian(nexcit,nexcit),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating hessian'
        endif
        allocate(gradient(nexcit),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating gradient'
        endif
        allocate(rrscr(lrscr),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating rrscr'
        endif
        allocate(hessian_evals(nexcit),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating rrscr'
        endif
        allocate(hessian_eigenvec(nexcit,nexcit),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating hessian_eigenvec'
        endif
        allocate(inversehessian(nexcit,nexcit),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating inversehessian'
        endif
        allocate(xfromh(nexcit),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating xfromh'
        endif
        allocate(test(nexcit,nexcit),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating test'
        endif

        ! Hessian iswritten in supermatrix form
        hessian(:,:) = 0.0_dp

        m = 0
        n = 0
        do l1=1,norb
            do l2=1,norb
                m = m + 1
                n = 0
                do l3=1,norb
                    do l4=1,norb
                        n = n + 1
                        hessian(n,m) = real(g2(l4,l3,l2,l1),dp)
                    enddo
                enddo
            enddo
        enddo
        
        if (m.ne.n) then
            stop 'Error: Hessian matrix is not in square form'
        endif

        if ((m.ne.nexcit).or.(n.ne.nexcit)) then
            stop 'Error: Hessian matrix has wrong number of elements'
        endif

        ! test hessian for hermiticity
        do l1=1,nexcit
            do l2=1,nexcit
                if (abs(hessian(l2,l1)-hessian(l1,l2)).gt.1e-8_dp) then
                    write(6,*) l2,l1,hessian(l2,l1),hessian(l1,l2)
                    stop 'Error: Hessian is not hermitian'
                endif
            enddo
        enddo

        ! gradient is written in supermatrix form
        gradient(:) = 0.0_dp

        m = 0
        do l1=1,norb
            do l2=1,norb
                m = m + 1
                gradient(m) = real(g1(l2,l1),dp)
            enddo
        enddo

        if ((m.ne.nexcit)) then
            stop 'Error: Gradient vector has wrong number of elements'
        endif


        write(6,*) 'The gradient vector is:'
        do l1=1,nexcit
            write(6,*) l1,gradient(l1)
        enddo

        write(6,*) 'The hessian vector is:'
        do l1=1,nexcit
            do l2=1,nexcit
                write(6,*) l2,l1,hessian(l2,l1)
            enddo
        enddo

        ! Hessian should be positive semi-definite 
        hessian_eigenvec = hessian
        call dsyev('V','U',nexcit,hessian_eigenvec,nexcit,&
            &hessian_evals,rrscr,lrscr,ierr)
        if (ierr.ne.0) then
            stop 'Error diagonalising Hessian'
        endif

        write(6,*) 'Eigenvectors of Hessian:'
        do l1=1,nexcit
            write(6,*) l1
            do l2=1,nexcit
                write(6,*) l2,hessian_eigenvec(l2,l1)
            enddo
        enddo

        write(6,*) 'Eigenvalues of Hessian:'
        do l1=1,nexcit
            write(6,*) l1,hessian_evals(l1)
        enddo

        mineval = abs(minval(hessian_evals)) + 2.0_dp


        ! x = (H-yI)^-1 g
        inversehessian(:,:) = 0.0_dp
        do l1=1,nexcit
            ! use (e + v)^-1 to deal with zero eigenvalues
            inversehessian(l1,l1) = 1.0_dp/(hessian_evals(l1) + mineval)
        enddo

        test(:,:) = 0.0_dp
        call dgemm('N','T',nexcit,nexcit,nexcit,1.0_dp,inversehessian,nexcit,&
        &hessian_eigenvec,nexcit,0.0_dp,test,nexcit)
        inversehessian(:,:) = 0.0_dp
        call dgemm('N','N',nexcit,nexcit,nexcit,1.0_dp,hessian_eigenvec,nexcit,&
        &test,nexcit,0.0_dp,inversehessian,nexcit)

        write(6,*) '----------------------------------------------------------'
        write(6,*) 'The inverse Hessian:'
        do l1=1,nexcit
            do l2=1,nexcit
                write(6,*) l2,l1,inversehessian(l2,l1)
            enddo
        enddo

        inversehessian = -1.0_dp*inversehessian

        xfromh(:) = 0.0_dp
        call dgemv('T',nexcit,nexcit,1.0_dp,inversehessian,nexcit,&
        &gradient,1,0.0_dp,xfromh,1)

       ! write(6,*) 'X-vector'
       ! do l1=1,nexcit
       !     if (abs(xfromh(l1)).gt.1e-12_dp) then
       !         write(6,*) l1,xfromh(l1)
       !     endif
       ! enddo


        cmat(:,:) = 0.0_dp
        m = 0
        do l1=1,norb
            do l2=1,norb
                m = m + 1
                cmat(l2,l1) = xfromh(m)
            enddo
        enddo

        write(6,*) '---------------------------------------------------------'
        write(6,*) 'The final transformation coefficients are:'
        do l1=1,norb
            do l2=1,norb
                write(6,*) l2,l1,cmat(l2,l1)
            enddo
        enddo

        write(6,*) '----------------------------------------------------------'

        deallocate(gradient)
        deallocate(hessian_eigenvec)
        deallocate(hessian_evals)
        deallocate(xfromh)
        deallocate(test)
        deallocate(inversehessian)
        deallocate(hessian)
        
    end subroutine NewtonRaphson

    subroutine PerformNewtonRaphson(cmat,g1,g2)

        ! This subroutine is to perform several Newton Raphson steps in a
        ! row in order to find a new minimum

        complex(dp), allocatable, intent(inout) :: cmat(:,:),g1(:,:)
        complex(dp), allocatable, intent(inout) :: g2(:,:,:,:)
        complex(dp), allocatable :: transcoeffs(:,:),originalcoeffs(:,:)
        complex(dp), allocatable :: xmat(:,:),dmat(:,:),e1(:,:)
        complex(dp), allocatable :: nos(:,:)
        real(dp), allocatable :: no_occ(:)
        real(dp), allocatable :: norrscr(:),nolrscr(:)
        complex(dp) :: energy,energy_old,sumgradient
        integer :: l1,l2,l3,l4,ierr,iter,rlscr

        rlscr = 4*norb
        allocate(transcoeffs(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating transcoeffs'
        endif
        allocate(originalcoeffs(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating originalcoeffs'
        endif
        allocate(xmat(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating xmat'
        endif
        allocate(dmat(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating dmat'
        endif
        allocate(e1(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating e1'
        endif
        allocate(norrscr(rlscr),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating norrscr'
        endif
        allocate(nolrscr(rlscr),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating nolrscr'
        endif
        allocate(nos(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating nos'
        endif
        allocate(no_occ(norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating no_occ'
        endif



        write(6,*) '------------------------------------------------------------'
        write(6,*) 'Performing a Newton-Raphson iteration scheme'
        write(6,*) '------------------------------------------------------------'


        originalcoeffs = cmat

        ! calculate original HF energy (assuming HF basis set)
        ! construct density matrix
        dmat(:,:) = 0.0_dp
        do l1=1,norb
            do l2=1,norb
                do l3=1,(nelec/2)
                    dmat(l2,l1) = dmat(l2,l1) + (2.0_dp*originalcoeffs(l1,l3)*&
                        &conjg(originalcoeffs(l2,l3)))
                enddo
            enddo
        enddo

        energy = 0.0_dp
        do l1=1,norb
            do l2=1,norb
                do l3=1,norb
                    do l4=1,norb
                        ! coulomb integral
                        energy = energy + (dmat(l4,l2)*dmat(l3,l1)*umat(l4,l3,l2,l1))
                        ! exchange integral
                        energy = energy - (0.5_dp*dmat(l4,l2)*&
                            &dmat(l3,l1)*umat(l4,l3,l1,l2))
                    enddo
                enddo
            enddo
        enddo
        energy = 0.5_dp*energy

        ! kinetic energy
        do l1=1,norb
            do l2=1,norb
                energy = energy + (dmat(l2,l1)*tmat(l2,l1))
            enddo
        enddo


        write(6,*) 'The original HF energy is:',energy
        energy_old = energy

        iter = 0
        transcoeffs = originalcoeffs

        do
          
            iter = iter + 1

            originalcoeffs = transcoeffs
            cmat = transcoeffs

            write(6,*) '--------------------------------------------------------------'
            write(6,*) 'Iteration:',iter
            write(6,*) 'Performing Newton-Raphson step...'
            ! Evaluate the gradient at this point
            call FindGradient(cmat,g1,g2)
            ! perform a Newton Raphson step
            call NewtonRaphson(cmat,g1,g2)
            xmat = cmat
            ! find rotation matrix
            call TransformFromGradient(cmat,xmat,-0.1_dp)

            ! perform rotation
            call zgemm('T','N',norb,norb,norb,(1.0_dp,0.0_dp),cmat,norb,&
            &originalcoeffs,norb,(0.0_dp,0.0_dp),transcoeffs,norb)

            write(6,*) 'New coefficients after the transformation'
            do l1=1,norb
                do l2=1,norb
                    write(6,*) l2,l1,transcoeffs(l2,l1)
                enddo
            enddo

    !        ! evaluate energy from E(X) = E(0) + (g^T X) + (1/2 X^T H X)
    !        e1(:,:) = 0.0_dp
    !        call zgemm('T','N',norb,norb,norb,(1.0_dp,0.0_dp),g1,norb,&
    !        &xmat,norb,(0.0_dp,0.0_dp),e1,norb)



            ! calculate enegy using closed shell expression
            ! E = \sum_i 2 <i|h|i> + \sum_ij ((2<ij|ij> - <ij|ji>)
            dmat(:,:) = 0.0_dp
            do l1=1,norb
                do l2=1,norb
                    do l3=1,(nelec/2)
                        dmat(l2,l1) = dmat(l2,l1) + (2.0_dp*transcoeffs(l1,l3)*&
                            &conjg(transcoeffs(l2,l3)))
                    enddo
                enddo
            enddo

            ! find natural orbitals
            nos = dmat
            call zheev('V','U',norb,nos,norb,no_occ,norrscr,rlscr,nolrscr,ierr)
            if (ierr.ne.0) then
                stop 'Error diagonalising nos'
            endif

            write(6,*) 'Natural Orbitals:'
            do l1=1,norb
                write(6,*) l1
                do l2=1,norb
                    write(6,*) l2,nos(l2,l1)
                enddo
            enddo

            write(6,*) 'Occupation numbers:'
            do l1=1,norb
                write(6,*) l1,no_occ(l1)
            enddo

            write(6,*) 'Sum of natural orbital occupation numbers:',sum(no_occ)

 
            energy = 0.0_dp
            do l1=1,norb
                do l2=1,norb
                    do l3=1,norb
                        do l4=1,norb
                            ! coulomb integral
                            energy = energy + (dmat(l4,l2)*dmat(l3,l1)*umat(l4,l3,l2,l1))
                            ! exchange integral
                            energy = energy - (0.5_dp*dmat(l4,l2)*&
                                &dmat(l3,l1)*umat(l4,l3,l1,l2))
                        enddo
                    enddo
                enddo
            enddo
            energy = 0.5_dp*energy

            ! kinetic energy
            do l1=1,norb
                do l2=1,norb
                    energy = energy + (dmat(l2,l1)*tmat(l2,l1))
                enddo
            enddo

            write(6,*) 'The HF energy of the previous iteration is:',energy_old
            write(6,*) 'The new HF energy is:', energy

            if (abs(energy_old-energy).lt.1e-8_dp) then
                write(6,*) 'HF energy of previous and current iteration:'
                write(6,*) energy_old,energy
                write(6,*) 'Convergence achieved !'
                write(6,*) 'Stopping calculation !'
                exit
            endif

            energy_old = energy

        enddo

        ! return coefficients
        cmat = transcoeffs

        deallocate(dmat)
        deallocate(e1)
        deallocate(xmat)
        deallocate(transcoeffs)
        deallocate(nos)
        deallocate(no_occ)
        deallocate(norrscr)
        deallocate(nolrscr)
        deallocate(originalcoeffs)

    end subroutine PerformNewtonRaphson

end program main
