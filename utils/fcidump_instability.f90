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
real(dp), allocatable :: scr(:),smat_eval(:),expsteps(:)
real(dp), allocatable :: expenergies(:)
complex(dp) :: expenergy
real(dp) :: ecore
complex(dp) :: cz
complex(dp), allocatable :: coeffs(:,:),temptrans(:,:),orighessian(:,:)
complex(dp), allocatable :: grade1(:,:),grade2(:,:,:,:),testcoeffs(:,:)
complex(dp), allocatable :: expdmat(:,:)
real(dp) :: z
integer :: loc(2),expminindex(1)
real(dp) :: ratio,angle,hf_energy,expmineval
real(dp), allocatable :: inst_all(:,:),transarray(:)
integer :: nmonoex,scheme,minimizer,followmode
integer(int64) :: orbsym(1000)
integer :: i,j,k,l,a,ierr,ispin,m,lscr
integer :: l1,l2,l3,l4,l5
real(dp), allocatable :: transform(:,:),fdiag(:)
integer, allocatable :: energyorder(:)
integer, allocatable :: monoexcitations(:,:)
logical :: exists,trealinds,tmolpro,complexint,error_diff
namelist /fci/ norb,nelec,ms2,orbsym,isym,iuhf,uhf,syml,symlz,propbitlen,nprop
    
    ! Defaults for rot_params input file
    trealinds = .true.
    tmolpro = .false.
    scheme = 3
    error_diff = .false.

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
    read(9,*) error_diff     ! this calculates the error in the gradient and hessian if they are evaluated exactly and using finite differences, the error is evaluated as a function of the size of the element in the rotation matrix (this helps to check that the gradient and hessian are correct)   
    read(9,*) minimizer      ! this option is to choose which minimzer is used in connection 
                             ! with scheme=3, 
                             ! if minimizer=1: an eigenvector-following scheme is applied 
                             ! to find a minimum, if minimizer=4: an eigenvector-following scheme
                             ! is applied to find a transition state
    read(9,*) followmode     ! this option is to choose which hessian mode should be followed if 
                             ! mode-following in an eigenvector-following scheme is applied in order 
                             ! to locate a transition state
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


    ! set up transformation matrix according to instability analysis

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
        allocate(orighessian(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating orighessian'
        endif
        allocate(expsteps(3200),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating expsteps'
        endif
        allocate(expenergies(3200),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating expenergies'
        endif
        allocate(testcoeffs(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating testcoeffs'
        endif
        allocate(expdmat(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating expdmat'
        endif


        ! this evaluates the gradient and hessian at the original point using the exact
        ! formulation and a finite difference method and evaluates the respective error 
        ! as a function of size of a particular element in the gradient/hessian matrix
        
        ! at the original HF solution
        if (error_diff) then
            coeffs(:,:) = 0.0_dp
            do l1=1,norb
                coeffs(l1,l1) = 1.0_dp
            enddo
            call ErrorGradientHessian(coeffs)
            call ErrorGradientHessian_2(coeffs)
        endif


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
                ratio = transarray(l1)
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
 
        ! this routine generates a starting guess wavefunction using a
        ! 'reverse perturbation' approach 
        !call StartingWavefunction(coeffs)

        ! these two routines would perform a standart RHF procedure either
        ! without or with DIIS
        ! the DIIS approach tends to lead back to the original solution
        !coeffs(:,:) = 0.0_dp
        !call RHF(coeffs,fdiag) 
        !call RHF_DIIS(coeffs,fdiag)

        
        ! this uses a Newton-Raphson iteration scheme in order to
        ! find a new minimum based on the previous instability analysis
        write(6,*) '----------------------------------------------------'
        write(6,*) 'Searching a new HF stationary point'
        write(6,*) 'based on instability analysis'
        write(6,*) '----------------------------------------------------'

        ! determine values of a/optimal step size for transformation
        ! where exp(-alpha U) where alpha is the step size
        write(6,*) '------------------------------------------------------'
        write(6,*) 'Searching the best value for the step size a'
        write(6,*) 'in the transformation exp(-aX)'
        write(6,*) '------------------------------------------------------'
        ! increase a in steps
        expsteps(:) = 0.0_dp
        expsteps(1) = -20.0_dp
        do l1=2,3200
            expsteps(l1) = expsteps(l1-1) + 0.01_dp
        enddo

        do l5=1,3200
            
            grade1 = coeffs
            temptrans(:,:) = 0.0_dp
            call TransformFromGradient(temptrans,grade1,expsteps(l5))
     
            grade1(:,:) = 0.0_dp
            do l1=1,norb
                grade1(l1,l1) = 1.0_dp
            enddo

            testcoeffs(:,:) = 0.0_dp
            call zgemm('T','N',norb,norb,norb,(1.0_dp,0.0_dp),temptrans,norb,&
            &grade1,norb,(0.0_dp,0.0_dp),testcoeffs,norb)

        
            ! calculate original HF energy (assuming HF basis set)
            ! construct density matrix
            expdmat(:,:) = 0.0_dp
            do l1=1,norb
                do l2=1,norb
                    do l3=1,(nelec/2)
                        expdmat(l2,l1) = expdmat(l2,l1) + (2.0_dp*testcoeffs(l1,l3)*&
                            &conjg(testcoeffs(l2,l3)))
                    enddo
                enddo
            enddo

            expenergy = 0.0_dp
            do l1=1,norb
                do l2=1,norb
                    do l3=1,norb
                        do l4=1,norb
                            ! coulomb integral
                            expenergy = expenergy + (expdmat(l4,l2)*expdmat(l3,l1)*umat(l4,l3,l2,l1))
                            ! exchange integral
                            expenergy = expenergy - (0.5_dp*expdmat(l4,l2)*&
                                        &expdmat(l3,l1)*umat(l4,l3,l1,l2))
                        enddo
                    enddo
                enddo
            enddo
            expenergy = 0.5_dp*expenergy

            ! kinetic energy
            do l1=1,norb
                do l2=1,norb
                    expenergy = expenergy + (expdmat(l2,l1)*tmat(l2,l1))
                enddo
            enddo

            expenergies(l5) = real(expenergy,dp) + ecore

        enddo

        ! chose a value
        if (minimizer.eq.2) then
            ! maximum for transition states
            expmineval = maxval(expenergies)
            expminindex = maxloc(expenergies)
        else
            ! minimum for other searches
            expmineval = minval(expenergies)
            expminindex = minloc(expenergies)
        endif


        writE(6,*) '---------------------------------------------------'
        write(6,*) 'Step sizes a and corresponding energies:'
        do l1=1,3200
            write(6,*) l1,expsteps(l1),expenergies(l1)
        enddo

        write(6,*) '-------------------------------------------------------'
        write(6,*) 'Choosing the following value for the step size:'
        write(6,*) (0+expminindex(1)),expsteps(0+expminindex(1)),expenergies(0+expminindex(1))
        write(6,*) '-------------------------------------------------------'
       
        ! perform rotation which gives a minimum energy
        grade1 = coeffs
        temptrans(:,:) = 0.0_dp
        call TransformFromGradient(temptrans,grade1,expsteps(0+expminindex(1)))
 
        grade1(:,:) = 0.0_dp
        do l1=1,norb
            grade1(l1,l1) = 1.0_dp
        enddo

        coeffs(:,:) = 0.0_dp
        call zgemm('T','N',norb,norb,norb,(1.0_dp,0.0_dp),temptrans,norb,&
        &grade1,norb,(0.0_dp,0.0_dp),coeffs,norb)

        grade2(:,:,:,:) = 0.0_dp

        if (minimizer.eq.1) then
            ! perform eigenvector-following
            ! to find a minimum
            call PerformNewtonRaphson_EigenVec_2(coeffs,grade1,grade2)
        elseif (minimizer.eq.2) then
            ! perform eigenvector-following
            ! to converge onto a transition state
            !call PerformNewtonRaphson_EigenVec_Trans(coeffs,grade1,grade2)
            call PerformNewtonRaphson_EigenVec_Trans_Mode(coeffs,grade1,grade2,followmode)
        endif

        ! as long as all coefficients are real this will be okay
        ! the orbitals should be real 
        transform = real(coeffs,dp)

        deallocate(coeffs)
        deallocate(temptrans)
        deallocate(grade1)
        deallocate(grade2)
        deallocate(orighessian)
        deallocate(expsteps)
        deallocate(testcoeffs)

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


        ! orbital energies
        ! only a rough guide since the fock matrix is not diagonal
        call DiagFockElements(fdiag)
        arr = fdiag


        ! new HF ground state energy approximately since
        ! fock matrix is no longer diagonal

        hf_energy = 0.0_dp
        if (uhf) then
            do l1=1,nelec
                hf_energy = hf_energy + tmat(l1,l1)
                do l2=1,nelec
                    hf_energy = hf_energy + umat(l1,l2,l1,l2)
                    if (mod(l1,2).eq.mod(l2,2)) then
                        hf_energy = hf_energy - umat(l1,l2,l2,l1)
                    endif
                enddo
            enddo
        elseif (.not.uhf) then
            do l1=1,(nelec/2)
                hf_energy = hf_energy + tmat(l1,l1)
                do l2=1,(nelec/2)
                    hf_energy = hf_energy + (2.0_dp*umat(l1,l2,l1,l2))&
                        & - umat(l1,l2,l2,l1)
                enddo
            enddo
        endif

        hf_energy = hf_energy + ecore

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
        real(dp), allocatable :: rrscr(:)
        complex(dp), allocatable :: rscr(:)
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

            energy = energy + ecore

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
            !if (iter.eq.1) then
            if ((abs(energy-energy_old).lt.1e-12_dp).and.&
                &(abs(std).lt.1e-12_dp)) then
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

            energy = energy + ecore

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
        real(dp), allocatable :: rrscr(:)
        complex(dp), allocatable :: rscr(:)
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

            energy = energy + ecore

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
            if ((abs(energy-energy_old).lt.1e-12_dp).and.&
                &(abs(std).lt.1e-12_dp)) then
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
        integer :: ierr

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
        complex(dp), allocatable :: grade2_test(:,:,:,:),grade2_test2(:,:,:,:)
        complex(dp), allocatable :: inactfmat(:,:),genfmatinactive(:,:)
        complex(dp), allocatable :: genfmat(:,:),ymat(:,:,:,:)
        complex(dp), allocatable, intent(inout) :: grade1(:,:),grade2(:,:,:,:)
        complex(dp), allocatable :: diffgrade21(:,:,:,:),diffgrade22(:,:,:,:)
        complex(dp) :: trace_1rdm,trace_2rdm
        complex(dp) :: energy_test
        integer :: l1,l2,l3,l4,l5,l6,ierr

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
        allocate(grade2_test(norb,norb,norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating grade2_test'
        endif
        allocate(grade2_test2(norb,norb,norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating grade2_test2'
        endif
        allocate(inactfmat(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating inactfmat'
        endif
        allocate(diffgrade21(norb,norb,norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating diffgrade21'
        endif
        allocate(diffgrade22(norb,norb,norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating diffgrade22'
        endif
        allocate(genfmatinactive(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating genfmatinactive'
        endif
 

        write(6,*) '----------------------------------------------------'
        write(6,*) 'Calculating electronic gradient and hessian'
        write(6,*) '-----------------------------------------------------'

        write(6,*) 'Constructing 1- and 2-RDM...'

        ! constructe 1-RDM for spatial orbitals. the sum over spin orbitals
        ! always gives a factor of 2 for RHF
        ! 1-RDM
        ! D_pq = <CSF|E_pq|CSF>
        d1mat(:,:) = 0.0_dp
        do l1=1,norb
            do l2=1,norb
                do l3=1,(nelec/2)
                    d1mat(l2,l1) = d1mat(l2,l1) + (2.0_dp*conjg(cmat(l2,l3))&
                        &*cmat(l1,l3))
                enddo
            enddo
        enddo

        ! The following is useful for debugging purposes
        ! for a single determinant closed shell system this way of calculating the 2-RDM
        ! should five the same result
        d2mat(:,:,:,:) = 0.0_dp
        do l1=1,norb
            do l2=1,norb
                do l3=1,norb
                    do l4=1,norb
                        d2mat(l4,l3,l2,l1) = ((d1mat(l4,l3)*d1mat(l2,l1)) -&
                            &(0.5_dp*d1mat(l4,l1)*d1mat(l2,l3)))
                    enddo
                enddo
            enddo
        enddo

        ! check: trace of 1-RDM should be nelec, trace 2-RDM should be 
        ! for 2-RDM in terms of spin orbitals
        ! nelec*(nelec-1) = (nelec^2-nelec) = \sum_pq d2mat(p,p,q,q) 
        ! - \sum_p d1mat(p,p)

        write(6,*) '-------------------------------------------------------'
        write(6,*) 'Diagonal elements of 1- and 2-RDM'
        write(6,*) '1-RDM: Occupation numbers'
        do l1=1,norb
            write(6,*) l1,d1mat(l1,l1)
        enddo
        write(6,*) '2-RDM: Simulateous occupation numbers'
        do l1=1,norb
            do l2=1,norb
                write(6,*) l2,l1,d2mat(l2,l2,l1,l1)
            enddo
        enddo
        write(6,*)

        trace_1rdm = 0.0_dp
        do l1=1,norb
            trace_1rdm = trace_1rdm + d1mat(l1,l1)
        enddo
        write(6,*) 'Trace of 1-RDM:',trace_1rdm
        
         
        trace_2rdm = 0.0_dp
        do l1=1,norb
            do l2=1,norb
                trace_2rdm = trace_2rdm + d2mat(l2,l2,l1,l1)
            enddo
        enddo
        write(6,*) 'Trace of 2-RDM:',trace_2rdm 


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

        ! calculate energy from general fock matrix in order to test it
        ! E(0) = (1/2)*\sum_pq (d_pq h_pq + delta_pq*f_pq)
        energy_test = 0.0_dp
        do l1=1,norb
            do l2=1,norb
                energy_test = energy_test + (d1mat(l2,l1)*tmat(l2,l1))
                if (l2.eq.l1) then
                    energy_test = energy_test + genfmat(l2,l1)
                endif
            enddo
        enddo

        energy_test = 0.5_dp*energy_test + ecore

        write(6,*) '---------------------------------------------------------'
        write(6,*) 'Energy obtained from generalised fock matrix:',energy_test

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
                        grade2(l4,l3,l2,l1) = grade2(l4,l3,l2,l1) - (2.0_dp*d1mat(l3,l2)&
                            &*tmat(l4,l1)) - (2.0_dp*ymat(l3,l4,l2,l1))
                        ! permute l2 and l1 
                        grade2(l4,l3,l2,l1) = grade2(l4,l3,l2,l1) - (2.0_dp*d1mat(l4,l1)&
                            &*tmat(l3,l2)) - (2.0_dp*ymat(l4,l3,l1,l2))
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
                                &+ (genfmat(l3,l2)+genfmat(l2,l3))
                        endif
                        if (l3.eq.l2) then
                            ! permute l2 and l1
                            grade2(l4,l3,l2,l1) = grade2(l4,l3,l2,l1) &
                                &+ (genfmat(l4,l1)+genfmat(l1,l4))
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
        deallocate(ymat)
        deallocate(grade2_test)
        deallocate(grade2_test2)
        deallocate(inactfmat)
        deallocate(diffgrade21)
        deallocate(diffgrade22)
        deallocate(genfmatinactive)

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
        integer :: l1,l2,ierr,rlscr


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
!
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
        integer :: l1,l2,ierr,rlscr


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


    subroutine NewtonRaphson_EigenVec_2(cmat,g1,g2)

        ! This subroutine is to apply an Eigenvector-Following approach in order to find a 
        ! minimum via a transformation x
        ! using the augmented Hessian and solving the eigenvalue problem
        ! (H   g) (X)           (x)
        ! (g^T 0) (1) = \lambda (1)
        ! where g is the gradient transformed into local hessian modes
        ! H the hessian in diagonal form
        ! such that the kth mode is being followed
        ! lambda_k = b_k/2 +/- 0.5*((b_k^2) + (4*f_k^2)^(1/2)
        ! where b_k is the eigenvalue of the kth mode which is being followed and 
        ! f_k is the respective element of the transformed gradient vector
        ! + is for maximisation and - for minimisation along the kth mode

        complex(dp), allocatable, intent(in) :: g1(:,:),g2(:,:,:,:)
        complex(dp), allocatable, intent(inout) :: cmat(:,:)
        real(dp), allocatable :: hessian(:,:),gradient(:)
        real(dp), allocatable :: rrscr(:),hessian_evals(:),hessian_eigenvec(:,:)
        real(dp), allocatable :: xfromh(:)
        real(dp), allocatable :: transformg(:,:)
        real(dp), allocatable :: aughessian(:,:),aughessian_eigenvec(:,:)
        real(dp), allocatable :: aughessian_evals(:),augrrscr(:)
        real(dp), allocatable :: inversehessian(:,:),temp(:)
        integer :: l1,l2,l3,l4,ierr,m,n,nexcit,lrscr,auglrscr

        ! number of matrix elements
        nexcit = norb*norb
        lrscr = 4*nexcit
        auglrscr = 4*(nexcit+1)

        write(6,*) '-----------------------------------------------------------'
        write(6,*) 'Eigenvector-Following Step for finding a minimum in the HF energy'
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
        allocate(aughessian((nexcit+1),(nexcit+1)),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating aughessian'
        endif
        allocate(aughessian_eigenvec((nexcit+1),(nexcit+1)),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating aughessian_eigenvec'
        endif
        allocate(aughessian_evals((nexcit+1)),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating aughessian_evals'
        endif
        allocate(augrrscr(auglrscr),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating augrrscr'
        endif
        allocate(transformg(nexcit,nexcit),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating transformg'
        endif
        allocate(temp(nexcit),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating temp'
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

        write(6,*) '--------------------------------------------------'

        ! Hessian should be positive semi-definite if the solution is to be stable
        ! diagonalise hessian
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


        temp(:) = 0.0_dp
        call dgemv('T',nexcit,nexcit,1.0_dp,hessian_eigenvec,nexcit,&
        &gradient,1,0.0_dp,temp,1)

        write(6,*) 'Gradient transformed into local hessian modes:'
        do l1=1,nexcit
            write(6,*) l1,temp(l1)
        enddo

        ! form augmented hessian
        aughessian(:,:) = 0.0_dp
        do l1=1,nexcit
            aughessian(l1,l1) = hessian_evals(l1)
            aughessian(l1,(nexcit+1)) = temp(l1)
            aughessian((nexcit+1),l1) = temp(l1)
        enddo

        aughessian_eigenvec = aughessian

        call dsyev('V','U',(nexcit+1),aughessian_eigenvec,(nexcit+1),&
            &aughessian_evals,augrrscr,auglrscr,ierr)
        if (ierr.ne.0) then
            stop 'Error diagonalising augmented Hessian'
        endif

        write(6,*) 'Eigenvectors of augmented Hessian:'
        do l1=1,(nexcit+1)
            write(6,'(I5)') l1
            do l2=1,(nexcit+1)
                write(6,*) l2,aughessian_eigenvec(l2,l1)
            enddo
        enddo

        write(6,*) 'Eigenvalues of augmented Hessian:'
        do l1=1,(nexcit+1)
            write(6,*) l1,aughessian_evals(l1)
        enddo

        ! calculate transformation vector
        ! x = - F_i V_i/(b_i - lambda)
        ! where V_i is the eigenvector of the hessian, b_i its eigenvalues
        ! lambda its lambda value and F_i the gradient component in local hessian 
        ! modes
        xfromh(:) = 0.0_dp
        do l1=1,nexcit
            do l2=1,nexcit
                if ((abs(hessian_eigenvec(l1,l2)).lt.1e-8_dp).or.&
                    &(abs(temp(l2)).lt.1e-8_dp)) cycle
                xfromh(l1) = xfromh(l1) + (temp(l2)*hessian_eigenvec(l1,l2)/&
                    &(aughessian_evals(1)-hessian_evals(l2)))
            enddo
        enddo


        ! return coefficients
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
        deallocate(rrscr)
        deallocate(hessian_eigenvec)
        deallocate(hessian_evals)
        deallocate(xfromh)
        deallocate(inversehessian)
        deallocate(aughessian)
        deallocate(aughessian_eigenvec)
        deallocate(aughessian_evals)
        deallocate(augrrscr)
        deallocate(temp)
        deallocate(hessian)
        
    end subroutine NewtonRaphson_EigenVec_2

    subroutine PerformNewtonRaphson_EigenVec_2(cmat,g1,g2)

        ! This subroutine is to perform several Newton Raphson steps in a
        ! row in order to find a new minimum
        ! this routine uses the eigenvector-following algorithm

        complex(dp), allocatable, intent(inout) :: cmat(:,:),g1(:,:)
        complex(dp), allocatable, intent(inout) :: g2(:,:,:,:)
        complex(dp), allocatable :: transcoeffs(:,:),originalcoeffs(:,:)
        complex(dp), allocatable :: xmat(:,:),dmat(:,:)!,e1(:,:)
        complex(dp), allocatable :: nos(:,:)
        real(dp), allocatable :: no_occ(:)
        complex(dp), allocatable :: norrscr(:)
        real(dp), allocatable :: nolrscr(:)
        real(dp) :: mineval
        real(dp) :: steps1(800),steps(800),stepsin(800),energies(800)
        integer :: minindex(1)
        complex(dp) :: energy,energy_old,energies1
        complex(dp) :: sumgradient,sumhessian
        integer :: l1,l2,l3,l4,l5,ierr,iter,rlscr

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
        write(6,*) 'Performing an Eigenvector-Following iteration scheme'
        write(6,*) 'for finding a minimum`'
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

        energy = energy + ecore

        write(6,*) 'The original HF energy is:',energy
        energy_old = energy

        ! if system has been moved away from HF solution
        steps1(:) = 0.0_dp
        steps1(1) = -5.0_dp
        do l1=2,800
            steps1(l1) = steps1(l1-1) + 0.1_dp
        enddo

        steps(:) = 0.0_dp
        steps(1) = -5.0_dp
        do l1=2,800
            steps(l1) = steps(l1-1) + 0.01_dp
        enddo


        !step = -0.1_dp
        iter = 0
        transcoeffs = originalcoeffs

        do
          
            iter = iter + 1

            originalcoeffs = transcoeffs
            cmat = transcoeffs

            if (iter.eq.1) then
                stepsin = steps1
            else
                stepsin = steps
            endif


            write(6,*) '--------------------------------------------------------------'
            write(6,*) 'Iteration:',iter
            write(6,*) 'Performing Eigenvector-Following step...'
            ! Evaluate the gradient at this point
            call FindGradient(cmat,g1,g2)
            ! calculate the sum of the gradient and hessian elements
            sumgradient = sum(abs(g1))
            sumhessian = sum(abs(g2))

            ! perform a Newton Raphson step
            call NewtonRaphson_EigenVec_2(cmat,g1,g2)
            !call NewtonRaphson(cmat,g1,g2)
            xmat = cmat

            do l5=1,800
                ! find rotation matrix
                call TransformFromGradient(cmat,xmat,stepsin(l5))
                
                ! perform rotation
                call zgemm('T','N',norb,norb,norb,(1.0_dp,0.0_dp),cmat,norb,&
                &originalcoeffs,norb,(0.0_dp,0.0_dp),transcoeffs,norb)

               ! calculate original HF energy (assuming HF basis set)
                ! construct density matrix
                dmat(:,:) = 0.0_dp
                do l1=1,norb
                    do l2=1,norb
                        do l3=1,(nelec/2)
                            dmat(l2,l1) = dmat(l2,l1) + (2.0_dp*transcoeffs(l1,l3)*&
                                &conjg(transcoeffs(l2,l3)))
                        enddo
                    enddo
                enddo

                energies1 = 0.0_dp
                do l1=1,norb
                    do l2=1,norb
                        do l3=1,norb
                            do l4=1,norb
                                ! coulomb integral
                                energies1 = energies1 + (dmat(l4,l2)*dmat(l3,l1)*umat(l4,l3,l2,l1))
                                ! exchange integral
                                energies1 = energies1 - (0.5_dp*dmat(l4,l2)*&
                                    &dmat(l3,l1)*umat(l4,l3,l1,l2))
                            enddo
                        enddo
                    enddo
                enddo
                energies1 = 0.5_dp*energies1

                ! kinetic energy
                do l1=1,norb
                    do l2=1,norb
                        energies1 = energies1 + (dmat(l2,l1)*tmat(l2,l1))
                    enddo
                enddo

                energies(l5) = real(energies1,dp) + ecore

            enddo

            ! chose a value
            mineval = minval(energies)
            minindex = minloc(energies)

            writE(6,*) '---------------------------------------------------'
            write(6,*) 'Step-sizes and energies'
            do l1=1,800
                write(6,*) l1,stepsin(l1),energies(l1)
            enddo

            write(6,*) 'Choosing the following value for the step-size:'
            write(6,*) minindex(1),stepsin(minindex(1)),energies(minindex(1))
            write(6,*) '-------------------------------------------------------'
            
            ! find rotation matrix
            call TransformFromGradient(cmat,xmat,stepsin(minindex(1)))

            ! perform rotation
            call zgemm('T','N',norb,norb,norb,(1.0_dp,0.0_dp),cmat,norb,&
            &originalcoeffs,norb,(0.0_dp,0.0_dp),transcoeffs,norb)

            write(6,*) 'New coefficients after the transformation:'
            do l1=1,norb
                do l2=1,norb
                    write(6,*) l2,l1,transcoeffs(l2,l1)
                enddo
            enddo


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

            write(6,*) '---------------------------------------------------------'
            write(6,*) 'Sum of natural orbital occupation numbers:',sum(no_occ)
            write(6,*) 'Maximum absolute element of gradient matrix:',&
                &maxval(abs(real(g1,dp)))
            write(6,*) 'Sum of absolute matrix elements of gradient and hessian:'
            write(6,*) sumgradient,sumhessian
            write(6,*) '----------------------------------------------------------'
 
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

            energy = energy + ecore

            write(6,*) 'The HF energy of the previous iteration:',energy_old
            write(6,*) 'The new HF energy:', energy


            if (abs(energy_old-energy).lt.1e-12_dp) then
                write(6,*) 'Convergence achieved !'
                write(6,*) 'Stopping calculation !'
                write(6,*) 'HF energy of previous and final iteration:'
                write(6,*) energy_old,energy
                write(6,*) 'Maximum residual absolute matrix element of gradient:',&
                    &maxval(abs(g1))
                write(6,*) 'Total sum of residual absolute value of gradient and hessian &
                    &matrix elements:'
                write(6,*) sumgradient,sumhessian
                exit
            endif

            energy_old = energy

        enddo

        ! return coefficients
        cmat = transcoeffs


        deallocate(dmat)
        deallocate(xmat)
        deallocate(transcoeffs)
        deallocate(nos)
        deallocate(no_occ)
        deallocate(norrscr)
        deallocate(nolrscr)
        deallocate(originalcoeffs)

    end subroutine PerformNewtonRaphson_EigenVec_2


    subroutine ErrorGradientHessian(coeffs)

        ! this subroutine calculates the error in the gradient when it is
        ! evaluated exactly and using a finite difference method, the error is evaluated
        ! as a function of the size of an element p,q in the gradient/hessian matrix

        complex(dp), allocatable, intent(in) :: coeffs(:,:)
        complex(dp), allocatable :: stepspq(:),temptrans(:,:)
        complex(dp), allocatable :: exactgrad(:,:),numgrad(:,:),numgradback(:,:)
        complex(dp), allocatable :: numgradcenter(:,:)
        complex(dp), allocatable :: gradcoeffs(:,:),gradcoeffsback(:,:)
        complex(dp), allocatable :: dmatorig(:,:),dmatnum(:,:),dmatnumback(:,:)
        complex(dp) :: energy_original,energy_num,energy_numback
        complex(dp), allocatable :: exacthess(:,:,:,:)
        complex(dp), allocatable :: energiesfor(:,:),energiesback(:,:)
        complex(dp), allocatable :: errorgradient(:,:),errorgradientback(:,:)
        complex(dp), allocatable :: errorgradientcenter(:,:)
        integer, allocatable :: p(:),q(:)
        integer :: l1,l2,l3,l4,l5,l6,l7,ierr,nsteps,m,n

        ! number of matrix elements to consider
        nsteps = (nelec/2)*(norb-(nelec/2))

        allocate(stepspq(235),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating stepspq'
        endif
        allocate(p(nelec/2),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating p'
        endif
        allocate(q(norb-(nelec/2)),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating q'
        endif
        allocate(exactgrad(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating exactgrad'
        endif
        allocate(numgrad(nsteps,235),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating numgrad'
        endif
        allocate(exacthess(norb,norb,norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating exacthess'
        endif
        allocate(temptrans(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating temptrans'
        endif
        allocate(gradcoeffs(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating gradcoeffs'
        endif
        allocate(dmatorig(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating dmatorig'
        endif
        allocate(dmatnum(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating dmatnum'
        endif
        allocate(errorgradient(nsteps,235),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating errorgradient'
        endif
        allocate(gradcoeffsback(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating gradcoeffsback'
        endif
        allocate(dmatnumback(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating dmatnum'
        endif
        allocate(numgradback(nsteps,235),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating numgradback'
        endif
        allocate(errorgradientback(nsteps,235),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating errorgradientback'
        endif
        allocate(energiesfor(nsteps,235),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating energiesfor'
        endif
        allocate(energiesback(nsteps,235),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating energiesback'
        endif
        allocate(numgradcenter(nsteps,235),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating numgradcenter'
        endif
        allocate(errorgradientcenter(nsteps,235),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating errorgradientcenter'
        endif
 

        write(6,*) '------------------------------------------------------------------'
        write(6,*) 'Evaluating error in gradient and hessian as a function of '
        write(6,*) 'size of matrix elements'
        write(6,*) '------------------------------------------------------------------'

        ! generate the steps by which the elements of the gradient/hessian matrix are 
        ! increased
        stepspq(:) = 1.0_dp
        stepspq(1) = 10.0_dp
        ! these points are only relevant for the behaviour of the energy
        do l1=2,91
            stepspq(l1) = stepspq(l1-1) - 0.1_dp
        enddo
        ! these points are relevant for the gradient
        do l1=92,235
            if (mod(l1,9).eq.2) then
                stepspq(l1) = 9.0_dp*stepspq(l1-1)/10.0_dp
            elseif (mod(l1,9).eq.3) then
                stepspq(l1) = 8.0_dp*stepspq(l1-1)/9.0_dp
            elseif (mod(l1,9).eq.4) then
                stepspq(l1) = 7.0_dp*stepspq(l1-1)/8.0_dp
            elseif (mod(l1,9).eq.5) then
                stepspq(l1) = 6.0_dp*stepspq(l1-1)/7.0_dp
            elseif (mod(l1,9).eq.6) then
                stepspq(l1) = 5.0_dp*stepspq(l1-1)/6.0_dp
            elseif (mod(l1,9).eq.7) then
                stepspq(l1) = 4.0_dp*stepspq(l1-1)/5.0_dp
            elseif (mod(l1,9).eq.8) then
                stepspq(l1) = 3.0_dp*stepspq(l1-1)/4.0_dp
            elseif (mod(l1,9).eq.0) then
                stepspq(l1) = 2.0_dp*stepspq(l1-1)/3.0_dp
            elseif (mod(l1,9).eq.1) then
                stepspq(l1) = stepspq(l1-1)/2.0_dp
            endif
        enddo

        write(6,*) 'Step sizes taken for each element'
        do l1=1,235
            write(6,*) l1,stepspq(l1)
        enddo

        ! set elements of gradien/hessian matrix whose step size should be incremented
        ! only non-redundant rotations need to be considered, i.e. those between occupied
        ! and unoccupied orbitals
        p(:) = 0
        q(:) = 0
        do l1=1,(nelec/2)
            p(l1) = l1
        enddo
        do l1=1,((norb-(nelec/2)))
            q(l1) = (nelec/2) + l1
        enddo

        ! calculate energy at this original point
        ! density matrix
        dmatorig(:,:) = 0.0_dp
        do l1=1,norb
            do l2=1,norb
                do l3=1,(nelec/2)
                    dmatorig(l2,l1) = dmatorig(l2,l1) + (2.0_dp*coeffs(l1,l3)*&
                        &conjg(coeffs(l2,l3)))
                enddo
            enddo
        enddo

        energy_original = 0.0_dp
        do l1=1,norb
            do l2=1,norb
                do l3=1,norb
                    do l4=1,norb
                        ! coulomb integral
                        energy_original = energy_original + (dmatorig(l4,l2)*dmatorig(l3,l1)*&
                            &umat(l4,l3,l2,l1))
                        ! exchange integral
                        energy_original = energy_original - (0.5_dp*dmatorig(l4,l2)*&
                            &dmatorig(l3,l1)*umat(l4,l3,l1,l2))
                    enddo
                enddo
            enddo
        enddo

        energy_original = 0.5_dp*energy_original + ecore

        ! kinetic energy
        do l1=1,norb
            do l2=1,norb
                energy_original = energy_original + (dmatorig(l2,l1)*tmat(l2,l1))
            enddo
        enddo

        write(6,*) '---------------------------------------------------------------'
        write(6,*) 'Original HF energy:',energy_original
        write(6,*) '----------------------------------------------------------------'

        write(6,*) 'Evaluating numerical gradient for an non-redunant rotations and &
            &step sizes'

        energiesfor(:,:) = 0.0_dp
        energiesback(:,:) = 0.0_dp

        m = 0
        do l5=1,(nelec/2)
            do l6=1,(norb-(nelec/2))
                m = m + 1
                do l7=1,235
                    ! steps in gradient
                    temptrans(:,:) = 0.0_dp
                    exactgrad(:,:) = 0.0_dp
                    exacthess(:,:,:,:) = 0.0_dp
                    exactgrad(p(l5),q(l6)) = stepspq(l7)
                    exactgrad(q(l6),p(l5)) = -stepspq(l7)

                    write(6,*) '----------------------------------------------------------'
                    write(6,*) 'Evaluating gradient for matrix elements',p(l5),q(l6)
                    write(6,*) 'using a step size of',stepspq(l7)

                    ! obtain transformation corresponding to gradient
                    call TransformFromGradient(temptrans,exactgrad,-1.0_dp)
                    gradcoeffs(:,:) = 0.0_dp
                    call zgemm('T','N',norb,norb,norb,(1.0_dp,0.0_dp),temptrans,norb,&
                        &coeffs,norb,(0.0_dp,0.0_dp),gradcoeffs,norb)

                    ! calculate energy at this new point
                    ! density matrix
                    dmatnum(:,:) = 0.0_dp
                    do l1=1,norb
                        do l2=1,norb
                            do l3=1,(nelec/2)
                                dmatnum(l2,l1) = dmatnum(l2,l1) + (2.0_dp*gradcoeffs(l1,l3)*&
                                    &conjg(gradcoeffs(l2,l3)))
                            enddo
                        enddo
                    enddo

                    energy_num = 0.0_dp
                    do l1=1,norb
                        do l2=1,norb
                            do l3=1,norb
                                do l4=1,norb
                                    ! coulomb integral
                                    energy_num = energy_num + (dmatnum(l4,l2)*dmatnum(l3,l1)*&
                                        &umat(l4,l3,l2,l1))
                                    ! exchange integral
                                    energy_num = energy_num - (0.5_dp*dmatnum(l4,l2)*&
                                        &dmatnum(l3,l1)*umat(l4,l3,l1,l2))
                                enddo
                            enddo
                        enddo
                    enddo

                    energy_num = 0.5_dp*energy_num + ecore

                    ! kinetic energy
                    do l1=1,norb
                        do l2=1,norb
                            energy_num = energy_num + (dmatnum(l2,l1)*tmat(l2,l1))
                        enddo
                    enddo

                    energiesfor(m,l7) = energy_num

                    write(6,*) 'Energies at forward point:',energy_num

                    ! the same for a backward transformation
                    ! obtain transformation corresponding to gradient
                    temptrans(:,:) = 0.0_dp
                    call TransformFromGradient(temptrans,exactgrad,1.0_dp)
                    gradcoeffsback(:,:) = 0.0_dp
                    call zgemm('T','N',norb,norb,norb,(1.0_dp,0.0_dp),temptrans,norb,&
                        &coeffs,norb,(0.0_dp,0.0_dp),gradcoeffsback,norb)

                    ! calculate energy at this new point
                    ! density matrix
                    dmatnumback(:,:) = 0.0_dp
                    do l1=1,norb
                        do l2=1,norb
                            do l3=1,(nelec/2)
                                dmatnumback(l2,l1) = dmatnumback(l2,l1) + (2.0_dp*gradcoeffsback(l1,l3)*&
                                    &conjg(gradcoeffsback(l2,l3)))
                            enddo
                        enddo
                    enddo

                    energy_numback = 0.0_dp
                    do l1=1,norb
                        do l2=1,norb
                            do l3=1,norb
                                do l4=1,norb
                                    ! coulomb integral
                                    energy_numback = energy_numback + (dmatnumback(l4,l2)*dmatnumback(l3,l1)*&
                                        &umat(l4,l3,l2,l1))
                                    ! exchange integral
                                    energy_numback = energy_numback - (0.5_dp*dmatnumback(l4,l2)*&
                                        &dmatnumback(l3,l1)*umat(l4,l3,l1,l2))
                                enddo
                            enddo
                        enddo
                    enddo

                    energy_numback = 0.5_dp*energy_numback + ecore

                    ! kinetic energy
                    do l1=1,norb
                        do l2=1,norb
                            energy_numback = energy_numback + (dmatnumback(l2,l1)*tmat(l2,l1))
                        enddo
                    enddo

                    energiesback(m,l7) = energy_numback

                    write(6,*) 'Energy backward point:',energy_numback

                    ! numerical gradient
                    ! using forward differencing
                    ! f'(x) = (f(x+h)-f(x))/h
                    ! ((E(0+deltak_pq)-E(0))/deltak_pq)
                    numgrad(m,l7) = (energy_num-energy_original)/stepspq(l7)

                    ! numerical gradient
                    ! using centered differencing
                    ! f'(x) = (f(x+h)-f(x-h))/(2h)
                    numgradcenter(m,l7) = (energy_num-energy_numback)/(2.0_dp*stepspq(l7))

                    !numerical gradient
                    ! using backward differencing
                    ! f'(x) = ((f(x)-f(x-h))/h
                    numgradback(m,l7) = (energy_original-energy_numback)/stepspq(l7)
               enddo
            enddo
        enddo

        write(6,*) '-----------------------------------------------------------------'
        write(6,*) 'The energies are: number,p,q,deltapq,E(0+deltapq),E(0-deltapq),E(0)'
        m = 0
        do l1=1,(nelec/2)
            do l2=1,(norb-(nelec/2))
                m = m + 1
                do l3=1,235
                    write(6,*) m,p(l1),q(l2),stepspq(l3),energiesfor(m,l3),&
                        &energiesback(m,l3),energy_original
                enddo
            enddo
        enddo



        write(6,*) '---------------------------------------------------------------'
        write(6,*) 'The numerical gradients are: number,p,q,deltapq,gradient_pq (forward, &
            &backward and centered differencing'
        m = 0
        do l1=1,(nelec/2)
            do l2=1,(norb-(nelec/2))
                m = m + 1
                do l3=1,235
                    write(6,*) m,p(l1),q(l2),stepspq(l3),numgrad(m,l3),numgradback(m,l3),&
                        &numgradcenter(m,l3)
                enddo
            enddo
        enddo

        ! calculating the exact gradient at the original solution
        call FindGradient(coeffs,exactgrad,exacthess)

        write(6,*) '------------------------------------------------------------------'
        write(6,*) 'The analytical gradient is:'
        m = 0
        do l1=1,(nelec/2)
            do l2=1,(norb-(nelec/2))
                m = m + 1
                write(6,*) p(l1),q(l2),exactgrad(p(l1),q(l2))
            enddo
        enddo
        

        m = 0
        do l1=1,(nelec/2)
            do l2=1,(norb-(nelec/2))
                m = m + 1
                do l3=1,235
                    errorgradient(m,l3) = numgrad(m,l3) - exactgrad(p(l1),q(l2))
                    errorgradientback(m,l3) = numgradback(m,l3) - exactgrad(p(l1),q(l2))
                    errorgradientcenter(m,l3) = numgradcenter(m,l3) -exactgrad(p(l1),q(l2))
                enddo
            enddo
        enddo

        write(6,*) '------------------------------------------------------------'
        write(6,*) 'The errors in the numerical gradient are: number,p,q,delta_pq,error &
            &(forward, backeard and centered differencing)'
        m = 0
        do l1=1,(nelec/2)
            do l2=1,(norb-(nelec/2))
                m = m + 1
                do l3=1,235
                    write(6,*) m,p(l1),q(l2),stepspq(l3),errorgradient(m,l3),&
                        &errorgradientback(m,l3),errorgradientcenter(m,l3)
                enddo
            enddo
        enddo

        ! write this out in a separate file for firther analysis
        open(9,file='Gradients_Error',status='unknown')
        write(9,*) '# Element p, element q, Increment delta_pq, Numerical gradient (forward,&
            & backward and centered differencing), Analytical&
            & gradient, Error in numerical gradient (forward, backward and centered &
            &differencing), Absolute error in numerical gradient (forward, backward&
            & and centered differencing)'
        m = 0
        do l1=1,(nelec/2)
            do l2=1,(norb-(nelec/2))
                m = m + 1
                do l3=1,235
                    write(9,'(2i5,19g20.12)') p(l1),q(l2),stepspq(l3),numgrad(m,l3),&
                        &numgradback(m,l3),numgradcenter(m,l3),exactgrad(p(l1),q(l2)),&
                        &errorgradient(m,l3),errorgradientback(m,l3),errorgradientcenter(m,l3),&
                        &abs(errorgradient(m,l3)),abs(errorgradientback(m,l3)),&
                        &abs(errorgradientcenter(m,l3))
                enddo
                write(9,*)
            enddo
        enddo
        close(9)

        ! write energies into a file
        
        open(9,file='Energies_Gradient',status='unknown')
        write(9,*) 'The energies are: number,p,q,deltapq,E(0+deltapq),E(0-deltapq),E(0)'
        ! and for negative steps
        n = 0
        do l1=(nelec/2),1,-1
            do l2=(norb-(nelec/2)),1,-1
                n = n + 1
                do l3=1,235
                    write(9,'(3i5,8g20.12)') n,p(l1),q(l2),-stepspq(l3),energiesback(n,l3),&
                        &energiesfor(n,l3),energy_original
                enddo
                write(9,'(3i5,8g20.12)') 0,p(l1),q(l2),(0.0_dp,0.0_dp),energy_original,&
                   &energy_original,energy_original
                do l3=235,1,-1
                    write(9,'(3i5,8g20.12)') n,p(l1),q(l2),stepspq(l3),energiesfor(n,l3),&
                        &energiesback(n,l3),energy_original
                enddo
                write(9,*)
            enddo
        enddo
        close(9)

        deallocate(stepspq)
        deallocate(p)
        deallocate(q)
        deallocate(energiesback)
        deallocate(energiesfor)
        deallocate(exactgrad)
        deallocate(exacthess)
        deallocate(numgrad)
        deallocate(numgradback)
        deallocate(numgradcenter)
        deallocate(temptrans)
        deallocate(gradcoeffs)
        deallocate(gradcoeffsback)
        deallocate(dmatnum)
        deallocate(dmatnumback)
        deallocate(dmatorig)
        deallocate(errorgradient)
        deallocate(errorgradientback)
        deallocate(errorgradientcenter)

    end subroutine ErrorGradientHessian

    subroutine ErrorGradientHessian_2(coeffs)

        ! this subroutine calculates the error in the hessian when it is
        ! evaluated exactly and using a finite difference method, the error is evaluated
        ! as a function of the size of an element p,q in the gradient/hessian matrix

        complex(dp), allocatable, intent(in) :: coeffs(:,:)
        complex(dp), allocatable :: stepspq(:),temptrans(:,:)
        complex(dp), allocatable :: exactgrad(:,:),exactgradin(:,:)
        complex(dp), allocatable :: exactgradstepback(:,:),exacthessstepback(:,:,:,:)
        complex(dp), allocatable :: exactgradstep(:,:),exacthessstep(:,:,:,:)
        complex(dp), allocatable :: numhess(:,:,:),numhessfor(:,:,:),numhessback(:,:,:)
        complex(dp), allocatable :: numhessgradback(:,:,:),numhessgrad(:,:,:)
        complex(dp), allocatable :: numhessgradcenter(:,:,:)
        complex(dp), allocatable :: gradcoeffsp1p1(:,:)
        complex(dp), allocatable :: dmatorig(:,:),dmatp1p1(:,:)
        complex(dp) :: energy_original,energy_p1p1,energy_p1m1,energy_m10,energy_0m1
        complex(dp) :: energy_m1p1,energy_m1m1,energy_0p1,energy_p10
        complex(dp) :: energy_p20,energy_m20
        complex(dp), allocatable :: exacthess(:,:,:,:)
        complex(dp), allocatable :: energiesp1p1(:,:,:),energiesp1m1(:,:,:)
        complex(dp), allocatable :: energiesm1p1(:,:,:),energiesm1m1(:,:,:)
        complex(dp), allocatable :: energiesp10(:,:,:),energies0p1(:,:,:)
        complex(dp), allocatable :: energiesm10(:,:,:),energies0m1(:,:,:)
        complex(dp), allocatable :: energiesp20(:,:,:),energiesm20(:,:,:)
        complex(dp), allocatable :: errorhess(:,:,:),errorhessgrad(:,:,:)
        complex(dp), allocatable :: errorhessfor(:,:,:),errorhessback(:,:,:)
        complex(dp), allocatable :: errorhessgradback(:,:,:)
        complex(dp), allocatable :: errorhessgradcenter(:,:,:)
        integer, allocatable :: p(:),q(:)
        integer :: l1,l2,l3,l4,l5,l6,l7,l8,l9,ierr,nsteps,m,n

        ! number of matrix elements to consider
        nsteps = (nelec/2)*(norb-(nelec/2))

        allocate(stepspq(235),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating stepspq'
        endif
        allocate(p(nelec/2),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating p'
        endif
        allocate(q(norb-(nelec/2)),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating q'
        endif
        allocate(exactgrad(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating exactgrad'
        endif
        allocate(numhess(nsteps,nsteps,235),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating numhess'
        endif
        allocate(numhessgrad(nsteps,nsteps,235),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating numhessgrad'
        endif
        allocate(numhessgradback(nsteps,nsteps,235),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating numhessgradback'
        endif
        allocate(numhessgradcenter(nsteps,nsteps,235),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating numhessgradcenter'
        endif
        allocate(exacthess(norb,norb,norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating exacthess'
        endif
        allocate(temptrans(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating temptrans'
        endif
        allocate(gradcoeffsp1p1(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating gradcoeffsp1p1'
        endif
        allocate(dmatorig(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating dmatorig'
        endif
        allocate(dmatp1p1(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating dmatp1p1'
        endif
        allocate(errorhessgrad(nsteps,nsteps,235),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating errorhessgrad'
        endif
        allocate(errorhess(nsteps,nsteps,235),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating errorhess'
        endif
        allocate(errorhessgradback(nsteps,nsteps,235),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating errorhessgradback'
        endif
        allocate(errorhessgradcenter(nsteps,nsteps,235),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating errorhessgradcenter'
        endif
        allocate(energiesp1p1(nsteps,nsteps,235),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating energiesp1p1'
        endif
        allocate(energiesp1m1(nsteps,nsteps,235),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating energiesp1m1'
        endif
        allocate(energiesm1p1(nsteps,nsteps,235),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating energiesm1p1'
        endif
        allocate(energiesm1m1(nsteps,nsteps,235),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating energiesm1m1'
        endif
        allocate(numhessback(nsteps,nsteps,235),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating numhessback'
        endif
        allocate(numhessfor(nsteps,nsteps,235),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating numhessfor'
        endif
        allocate(errorhessback(nsteps,nsteps,235),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating errorhessback'
        endif
        allocate(errorhessfor(nsteps,nsteps,235),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating errorhessfor'
        endif
        allocate(energiesp10(nsteps,nsteps,235),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating energiesp10'
        endif
        allocate(energies0p1(nsteps,nsteps,235),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating energies0p1'
        endif
        allocate(energiesm10(nsteps,nsteps,235),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating energiesm10'
        endif
        allocate(energies0m1(nsteps,nsteps,235),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating energies0m1'
        endif
        allocate(energiesp20(nsteps,nsteps,235),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating energiesp20'
        endif
        allocate(energiesm20(nsteps,nsteps,235),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating energiesm20'
        endif
        allocate(exactgradstep(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating exactgradstep'
        endif
         allocate(exacthessstep(norb,norb,norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating exacthessstep'
        endif
        allocate(exactgradstepback(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating exactgradstepback'
        endif
         allocate(exacthessstepback(norb,norb,norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating exacthessstepback'
        endif
        allocate(exactgradin(norb,norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating exactgradin'
        endif
 

        write(6,*) '------------------------------------------------------------------'
        write(6,*) 'Evaluating error in hessian as a function of '
        write(6,*) 'size of matrix elements'
        write(6,*) '------------------------------------------------------------------'

        ! generate the steps by which the elements of the gradient/hessian matrix are 
        ! increased
        stepspq(:) = 1.0_dp
        stepspq(1) = 10.0_dp
        ! these points are only relevant for the behaviour of the energy
        do l1=2,91
            stepspq(l1) = stepspq(l1-1) - 0.1_dp
        enddo
        ! these points are relevant for the hessian
        do l1=92,235
            if (mod(l1,9).eq.2) then
                stepspq(l1) = 9.0_dp*stepspq(l1-1)/10.0_dp
            elseif (mod(l1,9).eq.3) then
                stepspq(l1) = 8.0_dp*stepspq(l1-1)/9.0_dp
            elseif (mod(l1,9).eq.4) then
                stepspq(l1) = 7.0_dp*stepspq(l1-1)/8.0_dp
            elseif (mod(l1,9).eq.5) then
                stepspq(l1) = 6.0_dp*stepspq(l1-1)/7.0_dp
            elseif (mod(l1,9).eq.6) then
                stepspq(l1) = 5.0_dp*stepspq(l1-1)/6.0_dp
            elseif (mod(l1,9).eq.7) then
                stepspq(l1) = 4.0_dp*stepspq(l1-1)/5.0_dp
            elseif (mod(l1,9).eq.8) then
                stepspq(l1) = 3.0_dp*stepspq(l1-1)/4.0_dp
            elseif (mod(l1,9).eq.0) then
                stepspq(l1) = 2.0_dp*stepspq(l1-1)/3.0_dp
            elseif (mod(l1,9).eq.1) then
                stepspq(l1) = stepspq(l1-1)/2.0_dp
            endif
        enddo

        write(6,*) 'Step sizes taken for each element'
        do l1=1,235
            write(6,*) l1,stepspq(l1)
        enddo

        ! set elements of gradien/hessian matrix whose step size should be incremented
        ! only non-redundant rotations need to be considered, i.e. those between occupied
        ! and unoccupied orbitals
        p(:) = 0
        q(:) = 0
        do l1=1,(nelec/2)
            p(l1) = l1
        enddo
        do l1=1,((norb-(nelec/2)))
            q(l1) = (nelec/2) + l1
        enddo

        ! calculate energy at this original point
        ! density matrix
        dmatorig(:,:) = 0.0_dp
        do l1=1,norb
            do l2=1,norb
                do l3=1,(nelec/2)
                    dmatorig(l2,l1) = dmatorig(l2,l1) + (2.0_dp*coeffs(l1,l3)*&
                        &conjg(coeffs(l2,l3)))
                enddo
            enddo
        enddo

        energy_original = 0.0_dp
        do l1=1,norb
            do l2=1,norb
                do l3=1,norb
                    do l4=1,norb
                        ! coulomb integral
                        energy_original = energy_original + (dmatorig(l4,l2)*dmatorig(l3,l1)*&
                            &umat(l4,l3,l2,l1))
                        ! exchange integral
                        energy_original = energy_original - (0.5_dp*dmatorig(l4,l2)*&
                            &dmatorig(l3,l1)*umat(l4,l3,l1,l2))
                    enddo
                enddo
            enddo
        enddo

        energy_original = 0.5_dp*energy_original + ecore

        ! kinetic energy
        do l1=1,norb
            do l2=1,norb
                energy_original = energy_original + (dmatorig(l2,l1)*tmat(l2,l1))
            enddo
        enddo

        write(6,*) '---------------------------------------------------------------'
        write(6,*) 'Original HF energy:',energy_original
        write(6,*) '----------------------------------------------------------------'

        write(6,*) 'Evaluating numerical hessian for an non-redunant rotations and &
            &step sizes'

        energiesp1p1(:,:,:) = 0.0_dp
        energiesp1m1(:,:,:) = 0.0_dp
        energiesm1p1(:,:,:) = 0.0_dp
        energiesm1m1(:,:,:) = 0.0_dp
        energiesp10(:,:,:) = 0.0_dp
        energies0p1(:,:,:) = 0.0_dp
        energiesm10(:,:,:) = 0.0_dp
        energies0m1(:,:,:) = 0.0_dp

        m = 0
        do l5=1,(nelec/2)
            do l6=1,(norb-(nelec/2))
                m = m + 1
                n = 0
                do l8=1,(nelec/2)
                    do l9=1,(norb-(nelec/2))
                        n = n + 1
                        do l7=1,235
                            ! steps in gradient
                            temptrans(:,:) = 0.0_dp
                            exactgrad(:,:) = 0.0_dp
                            exacthess(:,:,:,:) = 0.0_dp
                            exactgrad(p(l5),q(l6)) = stepspq(l7)
                            exactgrad(q(l6),p(l5)) = -stepspq(l7)
                            exactgrad(p(l8),q(l9)) = stepspq(l7)
                            exactgrad(q(l9),p(l8)) = -stepspq(l7)

                            write(6,*) '----------------------------------------------------------'
                            write(6,*) 'Evaluating hessian for matrix elements',p(l5),q(l6),p(l8),q(l9)
                            write(6,*) 'using a step size of',stepspq(l7)

                            ! first at point +h1,+h2
                            ! obtain transformation corresponding to gradient
                            call TransformFromGradient(temptrans,exactgrad,-1.0_dp)
                            gradcoeffsp1p1(:,:) = 0.0_dp
                            call zgemm('T','N',norb,norb,norb,(1.0_dp,0.0_dp),temptrans,norb,&
                                &coeffs,norb,(0.0_dp,0.0_dp),gradcoeffsp1p1,norb)

                            ! calculate energy at this new point
                            ! density matrix
                            dmatp1p1(:,:) = 0.0_dp
                            do l1=1,norb
                                do l2=1,norb
                                    do l3=1,(nelec/2)
                                        dmatp1p1(l2,l1) = dmatp1p1(l2,l1) + (2.0_dp*gradcoeffsp1p1(l1,l3)*&
                                            &conjg(gradcoeffsp1p1(l2,l3)))
                                    enddo
                                enddo
                            enddo

                            energy_p1p1 = 0.0_dp
                            do l1=1,norb
                                do l2=1,norb
                                    do l3=1,norb
                                        do l4=1,norb
                                            ! coulomb integral
                                            energy_p1p1 = energy_p1p1 + (dmatp1p1(l4,l2)*dmatp1p1(l3,l1)*&
                                                &umat(l4,l3,l2,l1))
                                            ! exchange integral
                                            energy_p1p1 = energy_p1p1 - (0.5_dp*dmatp1p1(l4,l2)*&
                                                &dmatp1p1(l3,l1)*umat(l4,l3,l1,l2))
                                        enddo
                                    enddo
                                enddo
                            enddo

                            energy_p1p1 = 0.5_dp*energy_p1p1 + ecore

                            ! kinetic energy
                            do l1=1,norb
                                do l2=1,norb
                                    energy_p1p1 = energy_p1p1 + (dmatp1p1(l2,l1)*tmat(l2,l1))
                                enddo
                            enddo

                            energiesp1p1(m,n,l7) = energy_p1p1

                            write(6,*) 'Energies at point (+h1,+h2):',energy_p1p1

                            ! the same for a backward transformation at point (-h1,-h2)
                            ! obtain transformation corresponding to gradient
                            temptrans(:,:) = 0.0_dp
                            call TransformFromGradient(temptrans,exactgrad,1.0_dp)
                            gradcoeffsp1p1(:,:) = 0.0_dp
                            call zgemm('T','N',norb,norb,norb,(1.0_dp,0.0_dp),temptrans,norb,&
                                &coeffs,norb,(0.0_dp,0.0_dp),gradcoeffsp1p1,norb)

                            ! calculate energy at this new point
                            ! density matrix
                            dmatp1p1(:,:) = 0.0_dp
                            do l1=1,norb
                                do l2=1,norb
                                    do l3=1,(nelec/2)
                                        dmatp1p1(l2,l1) = dmatp1p1(l2,l1) + (2.0_dp*gradcoeffsp1p1(l1,l3)*&
                                            &conjg(gradcoeffsp1p1(l2,l3)))
                                    enddo
                                enddo
                            enddo

                            energy_m1m1 = 0.0_dp
                            do l1=1,norb
                                do l2=1,norb
                                    do l3=1,norb
                                        do l4=1,norb
                                            ! coulomb integral
                                            energy_m1m1 = energy_m1m1 + (dmatp1p1(l4,l2)*dmatp1p1(l3,l1)*&
                                                &umat(l4,l3,l2,l1))
                                            ! exchange integral
                                            energy_m1m1 = energy_m1m1 - (0.5_dp*dmatp1p1(l4,l2)*&
                                                &dmatp1p1(l3,l1)*umat(l4,l3,l1,l2))
                                        enddo
                                    enddo
                                enddo
                            enddo

                            energy_m1m1 = 0.5_dp*energy_m1m1 + ecore

                            ! kinetic energy
                            do l1=1,norb
                                do l2=1,norb
                                    energy_m1m1 = energy_m1m1 + (dmatp1p1(l2,l1)*tmat(l2,l1))
                                enddo
                            enddo

                            energiesm1m1(m,n,l7) = energy_m1m1

                            write(6,*) 'Energy at point (-h1,-h2):',energy_m1m1

                            ! steps in gradient
                            temptrans(:,:) = 0.0_dp
                            exactgrad(:,:) = 0.0_dp
                            exacthess(:,:,:,:) = 0.0_dp
                            exactgrad(p(l5),q(l6)) = stepspq(l7)
                            exactgrad(q(l6),p(l5)) = -stepspq(l7)
                            exactgrad(p(l8),q(l9)) = -stepspq(l7)
                            exactgrad(q(l9),p(l8)) = stepspq(l7)

                            ! first at point +h1,-h2
                            ! obtain transformation corresponding to gradient
                            call TransformFromGradient(temptrans,exactgrad,-1.0_dp)
                            gradcoeffsp1p1(:,:) = 0.0_dp
                            call zgemm('T','N',norb,norb,norb,(1.0_dp,0.0_dp),temptrans,norb,&
                                &coeffs,norb,(0.0_dp,0.0_dp),gradcoeffsp1p1,norb)

                            ! calculate energy at this new point
                            ! density matrix
                            dmatp1p1(:,:) = 0.0_dp
                            do l1=1,norb
                                do l2=1,norb
                                    do l3=1,(nelec/2)
                                        dmatp1p1(l2,l1) = dmatp1p1(l2,l1) + (2.0_dp*gradcoeffsp1p1(l1,l3)*&
                                            &conjg(gradcoeffsp1p1(l2,l3)))
                                    enddo
                                enddo
                            enddo

                            energy_p1m1 = 0.0_dp
                            do l1=1,norb
                                do l2=1,norb
                                    do l3=1,norb
                                        do l4=1,norb
                                            ! coulomb integral
                                            energy_p1m1 = energy_p1m1 + (dmatp1p1(l4,l2)*dmatp1p1(l3,l1)*&
                                                &umat(l4,l3,l2,l1))
                                            ! exchange integral
                                            energy_p1m1 = energy_p1m1 - (0.5_dp*dmatp1p1(l4,l2)*&
                                                &dmatp1p1(l3,l1)*umat(l4,l3,l1,l2))
                                        enddo
                                    enddo
                                enddo
                            enddo

                            energy_p1m1 = 0.5_dp*energy_p1m1 + ecore

                            ! kinetic energy
                            do l1=1,norb
                                do l2=1,norb
                                    energy_p1m1 = energy_p1m1 + (dmatp1p1(l2,l1)*tmat(l2,l1))
                                enddo
                            enddo

                            energiesp1m1(m,n,l7) = energy_p1m1

                            write(6,*) 'Energies at point (+h1,-h2):',energy_p1m1

                            ! the same for a backward transformation at point (-h1,+h2)
                            ! obtain transformation corresponding to gradient
                            temptrans(:,:) = 0.0_dp
                            call TransformFromGradient(temptrans,exactgrad,1.0_dp)
                            gradcoeffsp1p1(:,:) = 0.0_dp
                            call zgemm('T','N',norb,norb,norb,(1.0_dp,0.0_dp),temptrans,norb,&
                                &coeffs,norb,(0.0_dp,0.0_dp),gradcoeffsp1p1,norb)

                            ! calculate energy at this new point
                            ! density matrix
                            dmatp1p1(:,:) = 0.0_dp
                            do l1=1,norb
                                do l2=1,norb
                                    do l3=1,(nelec/2)
                                        dmatp1p1(l2,l1) = dmatp1p1(l2,l1) + (2.0_dp*gradcoeffsp1p1(l1,l3)*&
                                            &conjg(gradcoeffsp1p1(l2,l3)))
                                    enddo
                                enddo
                            enddo

                            energy_m1p1 = 0.0_dp
                            do l1=1,norb
                                do l2=1,norb
                                    do l3=1,norb
                                        do l4=1,norb
                                            ! coulomb integral
                                            energy_m1p1 = energy_m1p1 + (dmatp1p1(l4,l2)*dmatp1p1(l3,l1)*&
                                                &umat(l4,l3,l2,l1))
                                            ! exchange integral
                                            energy_m1p1 = energy_m1p1 - (0.5_dp*dmatp1p1(l4,l2)*&
                                                &dmatp1p1(l3,l1)*umat(l4,l3,l1,l2))
                                        enddo
                                    enddo
                                enddo
                            enddo

                            energy_m1p1 = 0.5_dp*energy_m1p1 + ecore

                            ! kinetic energy
                            do l1=1,norb
                                do l2=1,norb
                                    energy_m1p1 = energy_m1p1 + (dmatp1p1(l2,l1)*tmat(l2,l1))
                                enddo
                            enddo

                            energiesm1p1(m,n,l7) = energy_m1p1

                            write(6,*) 'Energy at point (-h1,+h2):',energy_m1p1

                            ! steps in gradient
                            temptrans(:,:) = 0.0_dp
                            exactgrad(:,:) = 0.0_dp
                            exacthess(:,:,:,:) = 0.0_dp
                            exactgrad(p(l5),q(l6)) = stepspq(l7)
                            exactgrad(q(l6),p(l5)) = -stepspq(l7)

                            ! first at point +h1,0
                            ! obtain transformation corresponding to gradient
                            call TransformFromGradient(temptrans,exactgrad,-1.0_dp)
                            gradcoeffsp1p1(:,:) = 0.0_dp
                            call zgemm('T','N',norb,norb,norb,(1.0_dp,0.0_dp),temptrans,norb,&
                                &coeffs,norb,(0.0_dp,0.0_dp),gradcoeffsp1p1,norb)

                            ! calculate energy at this new point
                            ! density matrix
                            dmatp1p1(:,:) = 0.0_dp
                            do l1=1,norb
                                do l2=1,norb
                                    do l3=1,(nelec/2)
                                        dmatp1p1(l2,l1) = dmatp1p1(l2,l1) + (2.0_dp*gradcoeffsp1p1(l1,l3)*&
                                            &conjg(gradcoeffsp1p1(l2,l3)))
                                    enddo
                                enddo
                            enddo

                            energy_p10 = 0.0_dp
                            do l1=1,norb
                                do l2=1,norb
                                    do l3=1,norb
                                        do l4=1,norb
                                            ! coulomb integral
                                            energy_p10 = energy_p10 + (dmatp1p1(l4,l2)*dmatp1p1(l3,l1)*&
                                                &umat(l4,l3,l2,l1))
                                            ! exchange integral
                                            energy_p10 = energy_p10 - (0.5_dp*dmatp1p1(l4,l2)*&
                                                &dmatp1p1(l3,l1)*umat(l4,l3,l1,l2))
                                        enddo
                                    enddo
                                enddo
                            enddo

                            energy_p10 = 0.5_dp*energy_p10 + ecore

                            ! kinetic energy
                            do l1=1,norb
                                do l2=1,norb
                                    energy_p10 = energy_p10 + (dmatp1p1(l2,l1)*tmat(l2,l1))
                                enddo
                            enddo

                            energiesp10(m,n,l7) = energy_p10

                            write(6,*) 'Energies at point (+h1,0):',energy_p10

                            ! first at point -h1,0
                            ! obtain transformation corresponding to gradient
                            call TransformFromGradient(temptrans,exactgrad,1.0_dp)
                            gradcoeffsp1p1(:,:) = 0.0_dp
                            call zgemm('T','N',norb,norb,norb,(1.0_dp,0.0_dp),temptrans,norb,&
                                &coeffs,norb,(0.0_dp,0.0_dp),gradcoeffsp1p1,norb)

                            ! calculate energy at this new point
                            ! density matrix
                            dmatp1p1(:,:) = 0.0_dp
                            do l1=1,norb
                                do l2=1,norb
                                    do l3=1,(nelec/2)
                                        dmatp1p1(l2,l1) = dmatp1p1(l2,l1) + (2.0_dp*gradcoeffsp1p1(l1,l3)*&
                                            &conjg(gradcoeffsp1p1(l2,l3)))
                                    enddo
                                enddo
                            enddo

                            energy_m10 = 0.0_dp
                            do l1=1,norb
                                do l2=1,norb
                                    do l3=1,norb
                                        do l4=1,norb
                                            ! coulomb integral
                                            energy_m10 = energy_m10 + (dmatp1p1(l4,l2)*dmatp1p1(l3,l1)*&
                                                &umat(l4,l3,l2,l1))
                                            ! exchange integral
                                            energy_m10 = energy_m10 - (0.5_dp*dmatp1p1(l4,l2)*&
                                                &dmatp1p1(l3,l1)*umat(l4,l3,l1,l2))
                                        enddo
                                    enddo
                                enddo
                            enddo

                            energy_m10 = 0.5_dp*energy_m10 + ecore

                            ! kinetic energy
                            do l1=1,norb
                                do l2=1,norb
                                    energy_m10 = energy_m10 + (dmatp1p1(l2,l1)*tmat(l2,l1))
                                enddo
                            enddo

                            energiesm10(m,n,l7) = energy_m10

                            write(6,*) 'Energies at point (-h1,0):',energy_m10

                            ! first at point +2h1,0
                            ! obtain transformation corresponding to gradient
                            call TransformFromGradient(temptrans,exactgrad,-2.0_dp)
                            gradcoeffsp1p1(:,:) = 0.0_dp
                            call zgemm('T','N',norb,norb,norb,(1.0_dp,0.0_dp),temptrans,norb,&
                                &coeffs,norb,(0.0_dp,0.0_dp),gradcoeffsp1p1,norb)

                            ! calculate energy at this new point
                            ! density matrix
                            dmatp1p1(:,:) = 0.0_dp
                            do l1=1,norb
                                do l2=1,norb
                                    do l3=1,(nelec/2)
                                        dmatp1p1(l2,l1) = dmatp1p1(l2,l1) + (2.0_dp*gradcoeffsp1p1(l1,l3)*&
                                            &conjg(gradcoeffsp1p1(l2,l3)))
                                    enddo
                                enddo
                            enddo

                            energy_p20 = 0.0_dp
                            do l1=1,norb
                                do l2=1,norb
                                    do l3=1,norb
                                        do l4=1,norb
                                            ! coulomb integral
                                            energy_p20 = energy_p20 + (dmatp1p1(l4,l2)*dmatp1p1(l3,l1)*&
                                                &umat(l4,l3,l2,l1))
                                            ! exchange integral
                                            energy_p20 = energy_p20 - (0.5_dp*dmatp1p1(l4,l2)*&
                                                &dmatp1p1(l3,l1)*umat(l4,l3,l1,l2))
                                        enddo
                                    enddo
                                enddo
                            enddo

                            energy_p20 = 0.5_dp*energy_p20 + ecore

                            ! kinetic energy
                            do l1=1,norb
                                do l2=1,norb
                                    energy_p20 = energy_p20 + (dmatp1p1(l2,l1)*tmat(l2,l1))
                                enddo
                            enddo

                            energiesp20(m,n,l7) = energy_p20

                            write(6,*) 'Energies at point (+2h1,0):',energy_p20

                            ! first at point -2h1,0
                            ! obtain transformation corresponding to gradient
                            call TransformFromGradient(temptrans,exactgrad,2.0_dp)
                            gradcoeffsp1p1(:,:) = 0.0_dp
                            call zgemm('T','N',norb,norb,norb,(1.0_dp,0.0_dp),temptrans,norb,&
                                &coeffs,norb,(0.0_dp,0.0_dp),gradcoeffsp1p1,norb)

                            ! calculate energy at this new point
                            ! density matrix
                            dmatp1p1(:,:) = 0.0_dp
                            do l1=1,norb
                                do l2=1,norb
                                    do l3=1,(nelec/2)
                                        dmatp1p1(l2,l1) = dmatp1p1(l2,l1) + (2.0_dp*gradcoeffsp1p1(l1,l3)*&
                                            &conjg(gradcoeffsp1p1(l2,l3)))
                                    enddo
                                enddo
                            enddo

                            energy_m20 = 0.0_dp
                            do l1=1,norb
                                do l2=1,norb
                                    do l3=1,norb
                                        do l4=1,norb
                                            ! coulomb integral
                                            energy_m20 = energy_m20 + (dmatp1p1(l4,l2)*dmatp1p1(l3,l1)*&
                                                &umat(l4,l3,l2,l1))
                                            ! exchange integral
                                            energy_m20 = energy_m20 - (0.5_dp*dmatp1p1(l4,l2)*&
                                                &dmatp1p1(l3,l1)*umat(l4,l3,l1,l2))
                                        enddo
                                    enddo
                                enddo
                            enddo

                            energy_m20 = 0.5_dp*energy_m20 + ecore

                            ! kinetic energy
                            do l1=1,norb
                                do l2=1,norb
                                    energy_m20 = energy_m20 + (dmatp1p1(l2,l1)*tmat(l2,l1))
                                enddo
                            enddo

                            energiesm20(m,n,l7) = energy_m20

                            write(6,*) 'Energies at point (-2h1,0):',energy_m20
 

                            ! steps in gradient
                            temptrans(:,:) = 0.0_dp
                            exactgrad(:,:) = 0.0_dp
                            exacthess(:,:,:,:) = 0.0_dp
                            exactgrad(p(l8),q(l9)) = stepspq(l7)
                            exactgrad(q(l9),p(l8)) = -stepspq(l7)

                            ! first at point 0,+h2
                            ! obtain transformation corresponding to gradient
                            call TransformFromGradient(temptrans,exactgrad,-1.0_dp)
                            gradcoeffsp1p1(:,:) = 0.0_dp
                            call zgemm('T','N',norb,norb,norb,(1.0_dp,0.0_dp),temptrans,norb,&
                                &coeffs,norb,(0.0_dp,0.0_dp),gradcoeffsp1p1,norb)

                            ! calculate energy at this new point
                            ! density matrix
                            dmatp1p1(:,:) = 0.0_dp
                            do l1=1,norb
                                do l2=1,norb
                                    do l3=1,(nelec/2)
                                        dmatp1p1(l2,l1) = dmatp1p1(l2,l1) + (2.0_dp*gradcoeffsp1p1(l1,l3)*&
                                            &conjg(gradcoeffsp1p1(l2,l3)))
                                    enddo
                                enddo
                            enddo

                            energy_0p1 = 0.0_dp
                            do l1=1,norb
                                do l2=1,norb
                                    do l3=1,norb
                                        do l4=1,norb
                                            ! coulomb integral
                                            energy_0p1 = energy_0p1 + (dmatp1p1(l4,l2)*dmatp1p1(l3,l1)*&
                                                &umat(l4,l3,l2,l1))
                                            ! exchange integral
                                            energy_0p1 = energy_0p1 - (0.5_dp*dmatp1p1(l4,l2)*&
                                                &dmatp1p1(l3,l1)*umat(l4,l3,l1,l2))
                                        enddo
                                    enddo
                                enddo
                            enddo

                            energy_0p1 = 0.5_dp*energy_0p1 + ecore

                            ! kinetic energy
                            do l1=1,norb
                                do l2=1,norb
                                    energy_0p1 = energy_0p1 + (dmatp1p1(l2,l1)*tmat(l2,l1))
                                enddo
                            enddo

                            energies0p1(m,n,l7) = energy_0p1

                            write(6,*) 'Energies at point (0,+h2):',energy_0p1

                            ! first at point 0,-h2
                            ! obtain transformation corresponding to gradient
                            call TransformFromGradient(temptrans,exactgrad,1.0_dp)
                            gradcoeffsp1p1(:,:) = 0.0_dp
                            call zgemm('T','N',norb,norb,norb,(1.0_dp,0.0_dp),temptrans,norb,&
                                &coeffs,norb,(0.0_dp,0.0_dp),gradcoeffsp1p1,norb)

                            ! calculate energy at this new point
                            ! density matrix
                            dmatp1p1(:,:) = 0.0_dp
                            do l1=1,norb
                                do l2=1,norb
                                    do l3=1,(nelec/2)
                                        dmatp1p1(l2,l1) = dmatp1p1(l2,l1) + (2.0_dp*gradcoeffsp1p1(l1,l3)*&
                                            &conjg(gradcoeffsp1p1(l2,l3)))
                                    enddo
                                enddo
                            enddo

                            energy_0m1 = 0.0_dp
                            do l1=1,norb
                                do l2=1,norb
                                    do l3=1,norb
                                        do l4=1,norb
                                            ! coulomb integral
                                            energy_0m1 = energy_0m1 + (dmatp1p1(l4,l2)*dmatp1p1(l3,l1)*&
                                                &umat(l4,l3,l2,l1))
                                            ! exchange integral
                                            energy_0m1 = energy_0m1 - (0.5_dp*dmatp1p1(l4,l2)*&
                                                &dmatp1p1(l3,l1)*umat(l4,l3,l1,l2))
                                        enddo
                                    enddo
                                enddo
                            enddo

                            energy_0m1 = 0.5_dp*energy_0m1 + ecore

                            ! kinetic energy
                            do l1=1,norb
                                do l2=1,norb
                                    energy_0m1 = energy_0m1 + (dmatp1p1(l2,l1)*tmat(l2,l1))
                                enddo
                            enddo

                            energies0m1(m,n,l7) = energy_0m1

                            write(6,*) 'Energies at point (0,-h2):',energy_0m1


                            ! numerical hessian
                            ! using centered differencing
                            ! d^2f/(dxdy) = (f(x+h1,y+h2)-f(x+h1,y-h2)-f(x-h1,y+h2)+f(x-h1,y-h2))/4h1h2
                            ! or if x=y
                            ! d^2f/(d^2x) = (f(x+h1,0)-2f(0,0)+f(x-h1,0))/(h1h2)
                            if ((p(l5).eq.p(l8)).and.(q(l6).eq.q(l9))) then
                                numhess(m,n,l7) = (energy_p10-(2.0_dp*energy_original)+energy_m10)&
                                    &/((stepspq(l7)**2))
                            else
 
                                numhess(m,n,l7) = (energy_p1p1-energy_p1m1-energy_m1p1+energy_m1m1)&
                                    &/(4.0_dp*(stepspq(l7)**2))
                            endif
 
                            ! numerical hessian
                            ! using forward differencing
                            ! d^2f/(dxdy) = (f(x+h1,y+h2)-f(x+h1,0)-f(0,y+h2)+f(0,0))/h1h2
                            ! or if x=y
                            ! d^2f/(d^2x) = (f(x+2h1,0)-2f(x+h1,0)+f(0,0))/(h1h2)
                            if ((p(l5).eq.p(l8)).and.(q(l6).eq.q(l9))) then
                                 numhessfor(m,n,l7) = (energy_p20-(2.0_dp*energy_p10)+energy_original)&
                                    &/((stepspq(l7)**2))
                            else
                                numhessfor(m,n,l7) = (energy_p1p1-energy_p10-energy_0p1+energy_original)&
                                    &/((stepspq(l7)**2))
                            endif
                             
                            ! numerical hessian
                            ! using backward differencing
                            ! d^2f/(dxdy) = (f(0,0)-f(x-h1,0)-f(0,y-h2)+f(x-h1,x-h2))/h1h2
                            ! or if x=y
                            ! d^2f/(d^2x) = (f(0,0)-2f(x-h1,0)+f(x-2h1,0))/(h1h2)
                            if ((p(l5).eq.p(l8)).and.(q(l6).eq.q(l9))) then
                                 numhessback(m,n,l7) = (energy_original-(2.0_dp*energy_m10)+energy_m20)&
                                    &/((stepspq(l7)**2))
                            else
                                numhessback(m,n,l7) = (energy_original-energy_m10-energy_0m1+energy_m1m1)&
                                    &/((stepspq(l7)**2))
                            endif
                        enddo
                    enddo
               enddo
            enddo
        enddo

        write(6,*) '-----------------------------------------------------------------'
        write(6,*) 'The energies are: number1,number2,p,q,r,s,deltapq/rs,E(0+deltapq,0+deltars),&
            &E(0+deltapq,0-deltars),E(0-deltapq,0+deltars),E(0-deltapq,0-deltars),&
            &E(0+deltapq,0),E(0,0+dealtapq),E(0-detapq,0),E(0,0-deltapq),E(0)'
        m = 0
        do l1=1,(nelec/2)
            do l2=1,(norb-(nelec/2))
                m = m + 1
                n = 0
                do l4=1,(nelec/2)
                    do l5=1,(norb-(nelec/2))
                        n = n + 1
                        do l3=1,235
                            write(6,*) m,n,p(l1),q(l2),p(l4),q(l5),stepspq(l3),energiesp1p1(m,n,l3),&
                                &energiesp1m1(m,n,l3),energiesm1p1(m,n,l3),&
                                &energiesm1m1(m,n,l3),energiesp10(m,n,l3),energies0p1(m,n,l3),&
                                &energiesm10(m,n,l3),energies0m1(m,n,l3),energy_original
                        enddo
                    enddo
                enddo
            enddo
        enddo


        write(6,*) '---------------------------------------------------------------'
        write(6,*) 'The numerical hessian is: number1,number2,p,q,r,s,deltapq/rs,&
            &hessian_pqrs (forward, backward and centered differencing)'
        m = 0
        do l1=1,(nelec/2)
            do l2=1,(norb-(nelec/2))
                m = m + 1
                n = 0
                do l4=1,(nelec/2)
                    do l5=1,(norb-(nelec/2))
                        n = n + 1
                        do l3=1,235
                            write(6,*) m,n,p(l1),q(l2),p(l4),q(l5),stepspq(l3),&
                                &numhessfor(m,n,l3),numhessback(m,n,l3),&
                                &numhess(m,n,l3)
                        enddo
                    enddo
                enddo
            enddo
        enddo

        ! calculating the exact and hessian at the original solution
        call FindGradient(coeffs,exactgrad,exacthess)

        write(6,*) '------------------------------------------------------------------'
        write(6,*) 'The analytical hessian is:'
        m = 0
        do l1=1,(nelec/2)
            do l2=1,(norb-(nelec/2))
                m = m + 1
                n = 0
                do l4=1,(nelec/2)
                    do l5=1,(norb-(nelec/2))
                        n = n + 1
                        write(6,*) p(l1),q(l2),p(l4),q(l5),exacthess(p(l1),q(l2),p(l4),q(l5))
                    enddo
                enddo
            enddo
        enddo
        

        m = 0
        do l1=1,(nelec/2)
            do l2=1,(norb-(nelec/2))
                m = m + 1
                n = 0
                do l4=1,(nelec/2)
                    do l5=1,(norb-(nelec/2))
                        n = n + 1
                        do l3=1,235
                            errorhess(m,n,l3) = numhess(m,n,l3) - exacthess(p(l1),q(l2),p(l4),q(l5))
                            errorhessfor(m,n,l3) = numhessfor(m,n,l3) - exacthess(p(l1),q(l2),p(l4),q(l5))
                            errorhessback(m,n,l3) = numhessback(m,n,l3) - exacthess(p(l1),q(l2),p(l4),q(l5))
                        enddo
                    enddo
                enddo
            enddo
        enddo

        write(6,*) '------------------------------------------------------------'
        write(6,*) 'The error in the numerical hessian is: number1,number2,,p,q,r,s,&
            &delta_pq/rs,error (forward, backward and centered differencing)'
        m = 0
        do l1=1,(nelec/2)
            do l2=1,(norb-(nelec/2))
                m = m + 1
                n = 0
                do l4=1,(nelec/2)
                    do l5=1,(norb-(nelec/2))
                        n = n + 1
                        do l3=1,235
                            write(6,*) m,n,p(l1),q(l2),p(l4),q(l5),stepspq(l3),&
                                &errorhessfor(m,n,l3),errorhessback(m,n,l3),&
                                &errorhess(m,n,l3)
                        enddo
                    enddo
                enddo
            enddo
        enddo

        write(6,*) '--------------------------------------------------------------------'
        write(6,*) 'Evaluating numerical hessian from gradient'
        write(6,*) '---------------------------------------------------------------------'

        ! find gradient at original solution
        call FindGradient(coeffs,exactgrad,exacthess)


        ! evaluate hessian numerically from the analytical gradient
        m = 0
        do l5=1,(nelec/2)
            do l6=1,(norb-(nelec/2))
                m = m + 1
                do l7=1,235
                    ! steps in gradient
                    temptrans(:,:) = 0.0_dp
                    exactgradin(:,:) = 0.0_dp
                    !exacthess(:,:,:,:) = 0.0_dp
                    exactgradin(p(l5),q(l6)) = stepspq(l7)
                    exactgradin(q(l6),p(l5)) = -stepspq(l7)

                    ! forward differencing
                    ! obtain transformation corresponding to gradient
                    call TransformFromGradient(temptrans,exactgradin,-1.0_dp)
                    gradcoeffsp1p1(:,:) = 0.0_dp
                    call zgemm('T','N',norb,norb,norb,(1.0_dp,0.0_dp),temptrans,norb,&
                    &coeffs,norb,(0.0_dp,0.0_dp),gradcoeffsp1p1,norb)

                    ! evaluate gradient at this step
                    call FindGradient(gradcoeffsp1p1,exactgradstep,exacthessstep)
 
                    ! backward differencing
                    ! obtain transformation corresponding to gradient
                    call TransformFromGradient(temptrans,exactgradin,1.0_dp)
                    gradcoeffsp1p1(:,:) = 0.0_dp
                    call zgemm('T','N',norb,norb,norb,(1.0_dp,0.0_dp),temptrans,norb,&
                    &coeffs,norb,(0.0_dp,0.0_dp),gradcoeffsp1p1,norb)

                    ! evaluate gradient at this step
                    call FindGradient(gradcoeffsp1p1,exactgradstepback,exacthessstepback)
                    
                    n = 0
                    do l8=1,(nelec/2)
                        do l9=1,(norb-(nelec/2))
                            n = n + 1
 
                            write(6,*) '----------------------------------------------------------'
                            write(6,*) 'Evaluating hessian from gradient for matrix elements',p(l5),q(l6),p(l8),q(l9)
                            write(6,*) 'using a step size of',stepspq(l7)

                            numhessgrad(m,n,l7) = (exactgradstep(p(l8),q(l9)) - &
                                &exactgrad(p(l8),q(l9)))/stepspq(l7)
                            errorhessgrad(m,n,l7) = numhessgrad(m,n,l7) - exacthess(p(l5),q(l6),p(l8),q(l9))
                            numhessgradback(m,n,l7) = (exactgrad(p(l8),q(l9)) - &
                                &exactgradstepback(p(l8),q(l9)))/stepspq(l7)
                            errorhessgradback(m,n,l7) = numhessgradback(m,n,l7) - exacthess(p(l5),q(l6),p(l8),q(l9))
                            numhessgradcenter(m,n,l7) = (exactgradstep(p(l8),q(l9)) - &
                                &exactgradstepback(p(l8),q(l9)))/(2.0_dp*stepspq(l7))
                            errorhessgradcenter(m,n,l7) = numhessgradcenter(m,n,l7) - exacthess(p(l5),q(l6),p(l8),q(l9))
                        enddo
                    enddo
                enddo
            enddo
        enddo

        write(6,*) '---------------------------------------------------------------'
        write(6,*) 'The numerical hessian from the gradient is: number1,number2,p,q,r,s,deltapq/rs,&
            &hessian_pqrs (forward, backward and centered differencing)'
        m = 0
        do l1=1,(nelec/2)
            do l2=1,(norb-(nelec/2))
                m = m + 1
                n = 0
                do l4=1,(nelec/2)
                    do l5=1,(norb-(nelec/2))
                        n = n + 1
                        do l3=1,235
                            write(6,*) m,n,p(l1),q(l2),p(l4),q(l5),stepspq(l3),&
                                &numhessgrad(m,n,l3),&
                                &numhessgradback(m,n,l3),numhessgradcenter(m,n,l3)
                        enddo
                    enddo
                enddo
            enddo
        enddo
        write(6,*) '------------------------------------------------------------'
        write(6,*) 'The error in the numerical hessian from the gradient is: &
            &number1,number2,,p,q,r,s,&
            &delta_pq/rs,error (forward, backward and centered differencing)'
        m = 0
        do l1=1,(nelec/2)
            do l2=1,(norb-(nelec/2))
                m = m + 1
                n = 0
                do l4=1,(nelec/2)
                    do l5=1,(norb-(nelec/2))
                        n = n + 1
                        do l3=1,235
                            write(6,*) m,n,p(l1),q(l2),p(l4),q(l5),stepspq(l3),&
                                &errorhessgrad(m,n,l3),&
                                &errorhessgradback(m,n,l3),errorhessgradcenter(m,n,l3)
                        enddo
                    enddo
                enddo
            enddo
        enddo

 
        ! write this out in a separate file for firther analysis
        open(9,file='Hessian_Error',status='unknown')
        write(9,*) '# Element p, Element q, Element r, Element s, Increment delta_pq/rs,&
            & Numerical hessian (forward, backward and centerd differencing), &
            & Numerical hessian from gradient (forward, backward and centered differencing), &
            & Analytical hessian,&
            & Error in numerical hessian (forward, backward and centered differencing),&
            & Error in numerical hessian from gradient (forward, backward and centered &
            &differencing), &
            & Absolute error in numerical hessian (forward, &
            &backward and centered differencing), Absolute error in numerical hessian from &
            &gradient (forward, backward and centered differencing)'
        m = 0
        do l1=1,(nelec/2)
            do l2=1,(norb-(nelec/2))
                m = m + 1
                n = 0
                do l4=1,(nelec/2)
                    do l5=1,(norb-(nelec/2))
                        n = n + 1
                        do l3=1,235
                            write(9,'(4i5,34g20.12)') p(l1),q(l2),p(l4),q(l5),stepspq(l3),&
                                &numhessfor(m,n,l3),numhessback(m,n,l3),numhess(m,n,l3),&
                                &numhessgrad(m,n,l3),numhessgradback(m,n,l3),&
                                &numhessgradcenter(m,n,l3),exacthess(p(l1),q(l2),p(l4),q(l5)),&
                                &errorhessfor(m,n,l3),errorhessback(m,n,l3),errorhess(m,n,l3),&
                                &errorhessgrad(m,n,l3),errorhessgradback(m,n,l3),&
                                &errorhessgradcenter(m,n,l3),abs(errorhessfor(m,n,l3)),&
                                &abs(errorhessback(m,n,l3)),abs(errorhess(m,n,l3)),&
                                &abs(errorhessgrad(m,n,l3)),abs(errorhessgradback(m,n,l3)),&
                                &abs(errorhessgradcenter(m,n,l3))
                        enddo
                        write(9,*)
                    enddo
                enddo
                write(9,*)
            enddo
        enddo
        close(9)

        deallocate(stepspq)
        deallocate(p)
        deallocate(q)
        deallocate(energiesp1p1)
        deallocate(energiesp1m1)
        deallocate(energiesm1p1)
        deallocate(energiesm1m1)
        deallocate(energies0p1)
        deallocate(energiesp10)
        deallocate(energiesm10)
        deallocate(energies0m1)
        deallocate(energiesp20)
        deallocate(energiesm20)
        deallocate(exactgradstepback)
        deallocate(exacthessstepback)
        deallocate(exactgrad)
        deallocate(exacthess)
        deallocate(exactgradstep)
        deallocate(exacthessstep)
        deallocate(exactgradin)
        deallocate(numhess)
        deallocate(numhessback)
        deallocate(numhessfor)
        deallocate(numhessgrad)
        deallocate(numhessgradback)
        deallocate(numhessgradcenter)
        deallocate(temptrans)
        deallocate(gradcoeffsp1p1)
        deallocate(dmatp1p1)
        deallocate(dmatorig)
        deallocate(errorhess)
        deallocate(errorhessback)
        deallocate(errorhessfor)
        deallocate(errorhessgrad)
        deallocate(errorhessgradback)
        deallocate(errorhessgradcenter)

    end subroutine ErrorGradientHessian_2


    subroutine NewtonRaphson_EigenVec_Trans_Mode(cmat,g1,g2,modein,origmode,switch)

        ! This subroutine is to apply an Eigenvector-Following approach in order to find a 
        ! transition state by maximizing along the kth mode via a transformation x
        ! using the augmented Hessian and solving the eigenvalue problem
        ! (H   g) (X)           (x)
        ! (g^T 0) (1) = \lambda (1)
        ! where g is the gradient transformed into local hessian modes
        ! H the hessian in diagonal form
        ! such that the kth mode is being followed
        ! lambda_k = b_k/2 +/- 0.5*((b_k^2) + (4*f_k^2)^(1/2)
        ! where b_k is the eigenvalue of the kth mode which is being followed and 
        ! f_k is the respective element of the transformed gradient vector
        ! + is for maximisation and - for minimisation along the kth mode

        complex(dp), allocatable, intent(in) :: g1(:,:),g2(:,:,:,:)
        complex(dp), allocatable, intent(inout) :: cmat(:,:)
        complex(dp), allocatable, intent(inout) :: origmode(:)
        logical, intent(in) :: switch
        integer :: mode,maxmodeloc(1)
        integer, intent(inout) :: modein
        complex(dp), allocatable :: overlap(:)
        real(dp), allocatable :: hessian(:,:),gradient(:)
        real(dp), allocatable :: rrscr(:),hessian_evals(:),hessian_eigenvec(:,:)
        real(dp), allocatable :: k_aughessian_evals(:),k_aughessian_eigenvec(:,:)
        real(dp), allocatable :: k_aughessian(:,:)
        real(dp), allocatable :: xfromh(:)
        real(dp), allocatable :: transformg(:,:)
        real(dp), allocatable :: aughessian(:,:),aughessian_eigenvec(:,:)
        real(dp), allocatable :: aughessian_evals(:),augrrscr(:)
        real(dp), allocatable :: inversehessian(:,:),temp(:)
        real(dp) :: maxmodeval
        integer :: l1,l2,l3,l4,ierr,m,n,nexcit,lrscr,auglrscr

        ! number of matrix elements
        nexcit = norb*norb
        lrscr = 4*nexcit
        auglrscr = 4*(nexcit+1)

        write(6,*) '-----------------------------------------------------------'
        write(6,*) 'Eigenvector-Following Step for finding a transition state &
            &in the HF energy'
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
            stop 'Error allocating hessian_evals'
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
        allocate(aughessian(nexcit,nexcit),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating aughessian'
        endif
        allocate(aughessian_eigenvec(nexcit,nexcit),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating aughessian_eigenvec'
        endif
        allocate(aughessian_evals((nexcit)),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating aughessian_evals'
        endif
        allocate(augrrscr(auglrscr),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating augrrscr'
        endif
        allocate(transformg(nexcit,nexcit),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating transformg'
        endif
        allocate(temp(nexcit),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating temp'
        endif
        allocate(k_aughessian_evals(2),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating k_aughessian_evals'
        endif
        allocate(k_aughessian_eigenvec(2,2),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating k_aughessian_eigenvec'
        endif
        allocate(k_aughessian(2,2),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating k_aughessian'
        endif
        allocate(overlap(nexcit),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating overlap'
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

        write(6,*) '--------------------------------------------------'

        ! Hessian should be positive semi-definite if the solution is to be stable
        ! diagonalise hessian
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


        temp(:) = 0.0_dp
        call dgemv('T',nexcit,nexcit,1.0_dp,hessian_eigenvec,nexcit,&
        &gradient,1,0.0_dp,temp,1)

        write(6,*) 'Gradient transformed into local hessian modes:'
        do l1=1,nexcit
            write(6,*) l1,temp(l1)
        enddo

        ! select the mode to be followed
        if (switch) then
            ! first iteration: choose one mode
            mode = modein
        else
            mode = modein
            overlap(:) = 0.0_dp
            do l1=1,nexcit
                overlap(l1) = dot_product(origmode,hessian_eigenvec(:,l1))
            enddo
            write(6,*) 'original folmode'
            do l1=1,nexcit
                write(6,*) l1,origmode(l1)
            enddo
            write(6,*) '--------------------------------------------------------'
            write(6,*) 'Overlap of present hessian modes with previously followed mode'
            do l1=1,nexcit
                write(6,*) l1,overlap(l1)
            enddo
            ! choose maximum overlap with mode of prevous cycle
            maxmodeval = maxval(abs(overlap(1:mode)))
            maxmodeloc = maxloc(abs(overlap(1:mode)))
            write(6,*) 'Choosing the following mode:',maxmodeloc(1),maxmodeval
            mode = maxmodeloc(1)
            modein = mode

        endif

        write(6,*) 'The followed mode is:',mode
        write(6,*) '-------------------------------------------------------'


        ! form augmented hessian
        aughessian(:,:) = 0.0_dp
        if (mode.eq.1) then
            do l1=2,nexcit
                aughessian((l1-1),(l1-1)) = hessian_evals(l1)
                aughessian((l1-1),nexcit) = temp(l1)
                aughessian(nexcit,(l1-1)) = temp(l1)
            enddo
        elseif (mode.eq.nexcit) then
            do l1=1,(nexcit-1)
                aughessian(l1,l1) = hessian_evals(l1)
                aughessian(l1,nexcit) = temp(l1)
                aughessian(nexcit,l1) = temp(l1)
            enddo
        else
            do l1=1,(mode-1)
                aughessian(l1,l1) = hessian_evals(l1)
                aughessian(l1,nexcit) = temp(l1)
                aughessian(nexcit,l1) = temp(l1)
            enddo
            do l1=(mode+1),nexcit
                aughessian((l1-1),(l1-1)) = hessian_evals(l1)
                aughessian((l1-1),nexcit) = temp(l1)
                aughessian(nexcit,(l1-1)) = temp(l1)
            enddo
        endif
 
        ! kth mode along which it is maximised
        ! choose lowest mode of hessian
        k_aughessian(:,:) = 0.0_dp
        k_aughessian(1,1) = hessian_evals(mode)
        k_aughessian(1,2) = temp(mode)
        k_aughessian(2,1) = temp(mode)

        aughessian_eigenvec = aughessian
        k_aughessian_eigenvec = k_aughessian

        call dsyev('V','U',nexcit,aughessian_eigenvec,nexcit,&
            &aughessian_evals,augrrscr,auglrscr,ierr)
        if (ierr.ne.0) then
            stop 'Error diagonalising augmented Hessian'
        endif

        call dsyev('V','U',2,k_aughessian_eigenvec,2,&
            &k_aughessian_evals,augrrscr,auglrscr,ierr)
        if (ierr.ne.0) then
            stop 'Error diagonalising augmented Hessian special mode'
        endif

        write(6,*) 'Eigenvectors of augmented Hessian:'
        do l1=1,nexcit
            write(6,'(I5)') l1
            do l2=1,nexcit
                write(6,*) l2,aughessian_eigenvec(l2,l1)
            enddo
        enddo

        write(6,*) 'Eigenvalues of augmented Hessian:'
        do l1=1,nexcit
            write(6,*) l1,aughessian_evals(l1)
        enddo

        write(6,*) 'Eigenvectors of augmented Hessian of special mode:'
        do l1=1,2
            write(6,'(I5)') l1
            do l2=1,2
                write(6,*) l2,k_aughessian_eigenvec(l2,l1)
            enddo
        enddo

        write(6,*) 'Eigenvalues of augmented Hessian of special mode:'
        do l1=1,2
            write(6,*) l1,k_aughessian_evals(l1)
        enddo

        ! calculate transformation vector
        ! x = - F_i V_i/(b_i - lambda)
        ! where V_i is the eigenvector of the hessian, b_i its eigenvalues
        ! lambda its lambda value and F_i the gradient component in local hessian 
        ! modes
        xfromh(:) = 0.0_dp
        do l1=1,nexcit
            if (mode.eq.1) then
                do l2=2,nexcit
                    if ((abs(hessian_eigenvec(l1,l2)).lt.1e-8_dp).or.&
                        &(abs(temp(l2)).lt.1e-8_dp)) cycle
                    xfromh(l1) = xfromh(l1) + (temp(l2)*hessian_eigenvec(l1,l2)/&
                        &(aughessian_evals(1)-hessian_evals(l2)))
                enddo
                if ((abs(hessian_eigenvec(l1,mode)).lt.1e-8_dp).or.&
                    &(abs(temp(mode)).lt.1e-1_dp)) cycle
                xfromh(l1) = xfromh(l1) - (temp(mode)*hessian_eigenvec(l1,mode)/&
                    &(k_aughessian_evals(2)-hessian_evals(mode)))
            elseif (mode.eq.nexcit) then
                do l2=1,(nexcit-1)
                    if ((abs(hessian_eigenvec(l1,l2)).lt.1e-8_dp).or.&
                        &(abs(temp(l2)).lt.1e-8_dp)) cycle
                    xfromh(l1) = xfromh(l1) + (temp(l2)*hessian_eigenvec(l1,l2)/&
                        &(aughessian_evals(1)-hessian_evals(l2)))
                enddo
                if ((abs(hessian_eigenvec(l1,mode)).lt.1e-8_dp).or.&
                    &(abs(temp(mode)).lt.1e-8_dp)) cycle
                xfromh(l1) = xfromh(l1) - (temp(mode)*hessian_eigenvec(l1,mode)/&
                    &(k_aughessian_evals(2)-hessian_evals(mode)))
            else
                do l2=1,(mode-1)
                    if ((abs(hessian_eigenvec(l1,l2)).lt.1e-8_dp).or.&
                        &(abs(temp(l2)).lt.1e-8_dp)) cycle
                    xfromh(l1) = xfromh(l1) + (temp(l2)*hessian_eigenvec(l1,l2)/&
                        &(aughessian_evals(1)-hessian_evals(l2)))
                enddo
                do l2=(mode+1),nexcit
                    if ((abs(hessian_eigenvec(l1,l2)).lt.1e-8_dp).or.&
                        &(abs(temp(l2)).lt.1e-8_dp)) cycle
                    xfromh(l1) = xfromh(l1) + (temp(l2)*hessian_eigenvec(l1,l2)/&
                        &(aughessian_evals(1)-hessian_evals(l2)))
                enddo
                if ((abs(hessian_eigenvec(l1,mode)).lt.1e-8_dp).or.&
                    &(abs(temp(mode)).lt.1e-8_dp)) cycle
                xfromh(l1) = xfromh(l1) - (temp(mode)*hessian_eigenvec(l1,mode)/&
                    &(k_aughessian_evals(2)-hessian_evals(mode)))
            endif
        enddo

        ! store the mode that has been followed
        origmode = hessian_eigenvec(:,mode)

        ! return coefficients
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
        deallocate(rrscr)
        deallocate(hessian_eigenvec)
        deallocate(hessian_evals)
        deallocate(k_aughessian)
        deallocate(k_aughessian_eigenvec)
        deallocate(k_aughessian_evals)
        deallocate(xfromh)
        deallocate(inversehessian)
        deallocate(aughessian)
        deallocate(aughessian_eigenvec)
        deallocate(aughessian_evals)
        deallocate(augrrscr)
        deallocate(temp)
        deallocate(hessian)
        deallocate(overlap)

    end subroutine NewtonRaphson_EigenVec_Trans_Mode

    subroutine PerformNewtonRaphson_EigenVec_Trans_Mode(cmat,g1,g2,followmodein)

        ! This subroutine is to perform several Newton Raphson steps in a
        ! row in order to find a new transition state
        ! this routine uses the eigenvector-following algorithm

        complex(dp), allocatable, intent(inout) :: cmat(:,:),g1(:,:)
        complex(dp), allocatable, intent(inout) :: g2(:,:,:,:)
        integer, intent(in) :: followmodein
        complex(dp), allocatable :: transcoeffs(:,:),originalcoeffs(:,:)
        complex(dp), allocatable :: xmat(:,:),dmat(:,:)
        complex(dp), allocatable :: nos(:,:)
        real(dp), allocatable :: no_occ(:)
        complex(dp), allocatable :: norrscr(:)
        real(dp), allocatable :: nolrscr(:)
        complex(dp), allocatable :: origmode(:)
        logical :: choosemode
        complex(dp) :: energy,energy_old
        complex(dp) :: sumgradient,sumhessian
        integer :: l1,l2,l3,l4,ierr,iter,rlscr,mode
        real(dp) :: sumgrad,sumgrad_old,step

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
        allocate(origmode(norb*norb),stat=ierr)
        if (ierr.ne.0) then
            stop 'Error allocating origmode'
        endif


        write(6,*) '------------------------------------------------------------'
        write(6,*) 'Performing an Eigenvector-Following iteration scheme'
        write(6,*) 'for finding a transition state'
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

        energy = energy + ecore

        write(6,*) 'The original HF energy is:',energy
        energy_old = energy

        !step = -0.1_dp
        iter = 0
        transcoeffs = originalcoeffs
        origmode = 0.0_dp

        ! mode to be followed 
        mode = followmodein 
        choosemode = .true.
        step = 0.1_dp

        do
          
            iter = iter + 1

            originalcoeffs = transcoeffs
            cmat = transcoeffs

            if (iter.eq.1) then
                choosemode = .true.
            else
                choosemode = .false.
            endif

            write(6,*) '--------------------------------------------------------------'
            write(6,*) 'Iteration:',iter
            write(6,*) 'Performing Eigenvector-Following step...'
            ! Evaluate the gradient at this point
            call FindGradient(cmat,g1,g2)
            ! calculate the sum of the gradient and hessian elements
            sumgradient = sum(abs(g1))
            sumhessian = sum(abs(g2))
            sumgrad = abs(sumgradient)
            if (iter.eq.1) sumgrad_old = abs(sumgradient)

            ! perform a Newton Raphson step
            call NewtonRaphson_EigenVec_Trans_Mode(cmat,g1,g2,mode,origmode,choosemode)
            !call NewtonRaphson(cmat,g1,g2)
            xmat = cmat

            ! find rotation matrix
            call TransformFromGradient(cmat,xmat,step)

            ! perform rotation
            call zgemm('T','N',norb,norb,norb,(1.0_dp,0.0_dp),cmat,norb,&
            &originalcoeffs,norb,(0.0_dp,0.0_dp),transcoeffs,norb)

            write(6,*) 'New coefficients after the transformation:'
            do l1=1,norb
                do l2=1,norb
                    write(6,*) l2,l1,transcoeffs(l2,l1)
                enddo
            enddo


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

            write(6,*) '---------------------------------------------------------'
            write(6,*) 'Sum of natural orbital occupation numbers:',sum(no_occ)
            write(6,*) 'Maximum absolute element of gradient matrix:',&
                &maxval(abs(real(g1,dp)))
            write(6,*) 'Sum of absolute matrix elements of gradient and hessian:'
            write(6,*) sumgradient,sumhessian
            write(6,*) '----------------------------------------------------------'
 
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

            energy = energy + ecore

            write(6,*) 'The HF energy of the previous iteration:',energy_old
            write(6,*) 'The new HF energy:', energy

            write(6,*) 'Sum of absolute gradient matrix elements (previous and current iteration):',sumgrad_old,sumgrad
            
            if ((abs(sumgrad).lt.0.5_dp)) then!
                write(6,*) 'Convergence achieved !'
                write(6,*) 'Stopping calculation !'
                write(6,*) 'HF energy of previous and final iteration:'
                write(6,*) energy_old,energy
                write(6,*) 'Maximum residual absolute matrix element of gradient:',&
                    &maxval(abs(g1))
                write(6,*) 'Total sum of residual absolute value of gradient and hessian &
                    &matrix elements:'
                write(6,*) sumgradient,sumhessian
                exit
            endif

            if (abs(sumgrad).lt.5.0_dp) then
                step = 0.01_dp
            endif

            energy_old = energy
            sumgrad_old = sumgrad

        enddo

        ! return coefficients
        cmat = transcoeffs


        deallocate(dmat)
        deallocate(xmat)
        deallocate(transcoeffs)
        deallocate(nos)
        deallocate(no_occ)
        deallocate(norrscr)
        deallocate(nolrscr)
        deallocate(originalcoeffs)
        deallocate(origmode)

    end subroutine PerformNewtonRaphson_EigenVec_Trans_Mode


end program main
