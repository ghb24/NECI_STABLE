! A new implementation file for csfs
module csf
    use systemdata, only: nel, brr, ecore, alat, nmsh, nbasismax, G1, nbasis
    use integralsdata, only: umat, fck, nmax
    use HElem
    use mt95, only: genrand_real2

    implicit none
    integer, parameter :: csf_orbital_mask = Z'1fffffff'
    integer, parameter :: csf_test_bit = 31
    integer, parameter :: csf_yama_bit = 30
    integer, parameter :: csf_ms_bit = 29

    ! Non-modularised functions (sigh)
    interface
        real*8 pure function choose(N,R)
            integer, intent(in) :: N,R
        end function
        logical function int_arr_eq (a, b, len)
            integer, intent(in), dimension(:) :: a, b
            integer, intent(in), optional :: len
        end function
        function GetHElement3_wrapper(NI,NJ,iC)
            use HElem
            use Systemdata, only: nEl
            INTEGER NI(nEl),NJ(nEl),iC
            type(HElement) GetHElement3_wrapper
        end function
        subroutine writedet(nunit,ni,nel,lterm)
            integer nunit,nel,ni(nel)
            logical lterm
        end subroutine
    end interface
    
contains
    ! Test if a specified determinant is a CSF
    logical pure function iscsf (nI)
        integer, dimension(:), intent(in) :: nI

        iscsf = btest(nI(1),csf_test_bit)
    end function

    ! This is (intentionally) a BRUTE FORCE way of calculating this.
    ! We would like to see if there is a nicer way of doing this using
    ! the representation matrices of the permutations which we are able 
    ! to calculate.
    type(HElement) function CSFGetHelement(NI, NJ)
        use memorymanager, only: LogMemAlloc, LogMemDealloc
        use SystemData, only: nel
        implicit none
        integer, intent(in) :: NI(nel), NJ(nel)
        integer nopen(2), nclosed(2), max_nopen, max_nup, max_ndets
        real*8 S(2), Ms(2)
        integer nup(2), ndets(2), i, j, det
        type(HElement) Hel, sum1

        integer, dimension (:), allocatable :: yama1, yama2
        real*8, dimension (:), allocatable :: coeffs1, coeffs2
        integer, dimension (:,:), allocatable :: dets1, dets2
        integer tagYamas(2)/0,0/, tagCoeffs(2)/0,0/, tagDets(2)/0,0/

        character(*), parameter :: this_routine = 'CSFGetHelement'

        call get_csf_data(NI, nel, nopen(1), nclosed(1), S(1), Ms(1))
        call get_csf_data(NJ, nel, nopen(2), nclosed(2), S(2), Ms(2))

        ! If S or Ms are not consistent, then return 0
        if ((S(1).ne.S(2)) .or. (Ms(1).ne.Ms(2))) then
            CSFGetHelement = HElement(0) 
            return
        endif

        ! Get electronic details
        ! Using S instead of Ms to calculate nup, as this has the fewest
        ! determinants, and the Ms=S case is degenerate.
        nup(1) = (nopen(1) + 2*S(1))/2
        ndets(1) = int(choose(nopen(1),nup(1)))
        nup(2) = (nopen(2) + 2*S(2))/2
        ndets(2) = int(choose(nopen(2),nup(2)))

        ! Allocate as required
        allocate (yama1(nopen(1)), yama2(nopen(2)))
        allocate (coeffs1(ndets(1)), coeffs2(ndets(2)))
        allocate (dets1(ndets(1), nel), dets2(ndets(2), nel))
        call LogMemAlloc ('yama1',size(yama1),4,this_routine,tagYamas(1))
        call LogMemAlloc ('yama2',size(yama2),4,this_routine,tagYamas(2))
        call LogMemAlloc ('coeffs1',size(coeffs1),8,this_routine,tagCoeffs(1))
        call LogMemAlloc ('coeffs2',size(coeffs2),8,this_routine,tagCoeffs(2))
        call LogMemAlloc ('dets1',size(dets1),4,this_routine,tagDets(1))
        call LogMemAlloc ('dets2',size(dets2),4,this_routine,tagDets(2))
        
        ! Calculate all possible permutations to construct determinants
        ! (Where 0=alpha, 1=beta when generating NI, NJ below)
        call csf_get_dets(nopen(1), nup(1), ndets(1), nel, dets1)
        if ((nopen(1).eq.nopen(2)) .and. (nup(1).eq.nup(2))) then
            dets2 = dets1
        else
            call csf_get_dets(nopen(2), nup(2), ndets(2), nel, dets2)
        endif

        ! Extract the Yamanouchi symbols from the CSFs
        call get_csf_yama (NI, yama1)
        call get_csf_yama (NJ, yama2)

        ! Get the coefficients
        do det=1,ndets(1)
            coeffs1(det) = csf_coeff(yama1,dets1(det,nclosed(1)+1:nel),&
                                     nopen(1))
        enddo
        do det=1,ndets(2)
            coeffs2(det) = csf_coeff(yama2,dets2(det,nclosed(2)+1:nel),&
                                     nopen(2))
        enddo

        ! Generate determinants from spatial orbitals specified in NI, NJ
        do det = 1,ndets(1)
            dets1(det,1:nclosed(1)) = iand(NI(1:nclosed(1)), csf_orbital_mask)
            dets1(det,nclosed(1)+1:nel) = &
                    csf_alpha_beta(NI(nclosed(1)+1:nel), &
                                   dets1(det,nclosed(1)+1:nel))
        enddo
        do det = 1,ndets(2)
            dets2(det,1:nclosed(2)) = iand(NJ(1:nclosed(2)), csf_orbital_mask)
            dets2(det,nclosed(2)+1:nel) = &
                    csf_alpha_beta(NJ(nclosed(2)+1:nel),&
                                   dets2(det,nclosed(2)+1:nel))
        enddo
        ! There will/may be faster ways of doing this
        call csf_sort_det_block (dets1, ndets(1), nopen(1))
        call csf_sort_det_block (dets1, ndets(2), nopen(2))

        ! TODO: implement symmetry if NI,NJ are the same except for yama
        CSFGetHelement = HElement(0)
        do i=1,ndets(1)
            sum1 = Helement(0)
            do j=1,ndets(2)
                Hel = GetHelement3_wrapper(dets1(i,:), dets2(j,:),-1)
                sum1 = sum1 + Hel * HElement(coeffs2(j))
            enddo
            CSFGetHelement = CSFGetHelement + sum1*HElement(coeffs1(i))
        enddo

        ! Deallocate for cleanup
        deallocate (coeffs1, coeffs2, yama1, yama2, dets1, dets2)
        call LogMemDealloc (this_routine, tagYamas(1))
        call LogMemDealloc (this_routine, tagYamas(2))
        call LogMemDealloc (this_routine, tagCoeffs(1))
        call LogMemDealloc (this_routine, tagCoeffs(2))
        call LogMemDealloc (this_routine, tagDets(1))
        call LogMemDealloc (this_routine, tagDets(2))
    end function



    ! Fill the last nopen electrons of each determinant with 0 (alpha) or
    ! 1 (beta) in all possible permutations with nup alpha electrons.
    subroutine csf_get_dets (nopen, nup, ndets, nel, dets)
        integer, intent(in) :: ndets, nup, nopen, nel
        integer, intent(out) :: dets (ndets, nel)
        integer comb(nup), i, j

        if (nopen.eq.0) return

        forall (i=1:nup) comb(i) = i
        dets(:,nel-nopen+1:) = 1
        do i=1,ndets
            forall (j=1:nup) dets(i,nel-nopen+comb(j)) = 0
            do j=1,nup
                if ((comb(j+1).ne.(comb(j)+1)) .or. (j.eq.nup)) then
                    comb(j) = comb(j) + 1
                    exit
                else
                    comb(j) = j
                endif
            enddo
        enddo
    end subroutine

    subroutine csf_sort_det_block (dets, ndets, nopen)
        integer, intent(in) :: ndets, nopen
        integer, intent(inout) :: dets (ndets,nel)
        integer :: i, open_pos, tmp_dets (ndets, nel), nclosed, npos

        tmp_dets = dets
        nclosed = nel - nopen
        npos = 1
        open_pos = nclosed + 1
        ! Closed e- listed in pairs --> can jump pairs
        do i=1,nclosed-1,2
            do while ( (open_pos .le. nel) .and. &
                       (tmp_dets(1,open_pos) .lt. tmp_dets(1,i)) )
                dets(:,npos) = tmp_dets(:,open_pos)
                open_pos = open_pos + 1
                npos = npos + 1
            enddo
            dets(:,npos:npos+1) = tmp_dets(:,i:i+1)
            npos = npos + 2
        enddo
        ! Any remaining open electrons will be untouched.            
    end subroutine

    ! For all of the open shell electrons, generate a list where
    ! 0=alpha, 1=beta for the specified determinant.
    function csf_get_dorder (NI, nel, nopen)
        integer csf_get_dorder(nopen)
        integer, intent(in) :: NI(nel)
        integer, intent(in) :: nopen, nel
        integer nclosed

        nclosed = nel - nopen
        csf_get_dorder = -(mod(NI(nclosed+1:nel),2)-1)
    end function


    subroutine csf_get_yamas (nopen, sfinal, yama, ncsf_max)
        real*8, intent(in) :: sfinal
        integer, intent(in) :: nopen, ncsf_max
        integer, intent(out) :: yama (ncsf_max, nopen)
        real*8 spin (ncsf_max, nopen)
        integer npos, csf, ncsf, ncsf_next

        if (nopen == 0) return

        spin(1,nopen) = sfinal
        ncsf = 1
        ncsf_next = ncsf
        do npos = nopen, 2, -1
            do csf=1,ncsf
                if (2*spin(csf,npos) .lt. npos) then
                    spin(csf,npos-1) = spin(csf,npos) + 0.5
                    yama(csf,npos) = 2
                    if (spin(csf,npos) .ne. 0) then
                        ncsf_next = ncsf_next + 1
                        if (ncsf_next .gt. ncsf_max) exit
                        !spin(ncsf_next,npos:nopen) = spin(csf,npos:nopen)
                        yama(ncsf_next,npos+1:nopen) = yama(csf,npos+1:nopen)
                        spin(ncsf_next,npos-1) = spin(csf,npos) - 0.5
                        yama(ncsf_next,npos) = 1
                    endif
                else
                    spin(csf,npos-1) = spin(csf,npos) - 0.5
                    yama(csf,npos) = 1
                endif
            enddo
            ncsf = ncsf_next
            if (ncsf .gt. ncsf_max) exit
        enddo
        yama(:,1) = 1
    end subroutine
        
    ! Convert num (a member of CI) to an alpha or beta spin orbital where
    ! det==0 --> alpha
    integer elemental function csf_alpha_beta (num, det)
        integer, intent(in) :: num, det

        csf_alpha_beta = iand(num, csf_orbital_mask) - 1
        if (det .eq. 1) then
            csf_alpha_beta = ior(csf_alpha_beta, 1)
        else
            csf_alpha_beta = iand(csf_alpha_beta,Z'fffffffe')
        endif
        csf_alpha_beta = csf_alpha_beta + 1
    end function

    ! Extract the Yamanouchi symbol from the supplied CSF
    ! This assumes that the passed yama array is the correct size, if
    ! it is too small, then the symbol will be truncated.
    subroutine get_csf_yama(NI, yama)
        integer, intent(in), dimension(:) :: NI
        integer, intent(out), dimension(:) :: yama
        integer i, nopen
        logical open_shell

        nopen = 0
        open_shell = .false.
        do i=1,size(NI)
            if ( (.not.open_shell) .and. btest(NI(i), csf_yama_bit)) then
                open_shell = .true.
            endif

            if (open_shell) then
                nopen = nopen + 1
                if (btest(NI(i), csf_yama_bit)) then
                    yama(nopen) = 1
                else
                    yama(nopen) = 2
                endif
                if (nopen == size(yama)) exit
            endif
        enddo
    end subroutine

    ! Obtains the number of open shell electrons, the total spin
    ! and the Ms value for the specified csf.
    subroutine get_csf_data(NI, nel, nopen, nclosed, S, Ms)
        integer, intent(in) :: NI(nel), nel
        integer, intent(out) :: nopen, nclosed
        real*8, intent(out) :: S, Ms
        integer i
        logical open_shell

        nopen = 0
        nclosed = 0
        S = 0
        Ms = 0
        open_shell = .false.
        do i=1,nel
            ! Closed shell until they are no longer paired.
            if ((.not.open_shell) .and. &
                btest(NI(i),csf_yama_bit)) open_shell = .true.
                
            if (.not. open_shell) then
                nclosed = nclosed + 1
            else
                if (btest(NI(i), csf_yama_bit)) then
                    S = S + 0.5
                else
                    S = S - 0.5
                endif
                if (btest(NI(i), csf_ms_bit)) then
                    Ms = Ms + 0.5
                else
                    Ms = Ms - 0.5
                endif
            endif
        enddo
        nopen = nel - nclosed
    end subroutine

    ! Calculates the total number of CSFs possible for a system with
    ! nOpen unpaired electrons, and a total spin of S.
    ! This is the same as the number of available Serber functions
    integer pure function get_num_csfs (nOpen, S)
        integer, intent(in) :: nOpen
        real*8, intent(in) :: S
        integer :: S2
        S2 = 2*S

        if ((nopen < 0) .or. (mod(nOpen+S2, 2) /= 0))then
            get_num_csfs = 0
        else
            get_num_csfs = (2*S2 + 2) * choose(nOpen, (nOpen+S2)/2)
            get_num_csfs = get_num_csfs / (nOpen + S2 + 2)
        endif
    end function

    ! TODO: This can be optimised (don't need to generate them all)
    ! TODO: Generate random by random branching perhaps (lots of genrands...)
    subroutine csf_apply_random_yama (nI, nopen, S, ncsf, tForceChange)
        integer, intent(inout) :: nI(nel)
        integer, intent(in) :: nopen
        integer, intent(out) :: ncsf
        real*8, intent(in) :: S
        logical, intent(in) :: tForceChange
        integer :: yamas (0:get_num_csfs(nopen, S), nopen), num
        real*8 :: r

        ! Generate the Yamanouchi Symbols
        ncsf = size(yamas(:,1))-1
        call csf_get_yamas (nopen, S, yamas(1:,:), ncsf)

        if (tForceChange .and. ncsf > 1) then
            call get_csf_yama (nI, yamas(0,:))
        endif

        ! Pick and apply a random one
        do while (.true.)
            call genrand_real2(r)
            num = int(r*ncsf) + 1
            if ((.not.tForceChange) .or. (ncsf<2) .or. &
                .not.int_arr_eq(yamas(num,:),yamas(0,:))) then

                call csf_apply_yama (nI, yamas(num, :))
                exit
            endif
        enddo
    end subroutine

    ! Apply a Yamanouchi symbol to a csf
    subroutine csf_apply_yama (NI, csf)
        integer, intent(in), dimension(:) :: csf
        integer, intent(inout) :: NI(nel)
        integer i

        NI = ibset(NI, csf_test_bit)
        do i=1,size(csf)
            if (csf(size(csf)-i+1) .eq. 1) then
                NI(nel-i+1) = ibset(NI(nel-i+1), csf_yama_bit)
            else
                NI(nel-i+1) = ibclr(NI(nel-i+1), csf_yama_bit)
            endif
        enddo
    end subroutine

    ! Apply a specified Ms value to the csf
    subroutine csf_apply_ms (NI, Ms, nopen)
        integer, intent(inout) :: NI(nel)
        integer, intent(in) :: nopen
        real*8, intent(in) :: Ms
        integer i, nup, ndown
        
        ndown = (nopen - 2*MS)/2
        do i=1,ndown
            NI(nel-i+1) = ibclr(NI(nel-i+1), csf_ms_bit)
        enddo
        do i=ndown+1,nopen
            NI(nel-i+1) = ibset(NI(nel-i+1), csf_ms_bit)
        enddo
    end subroutine

    ! Calculates the total spin from a csf
    real*8 function csf_spin (csf)
        integer, intent(in), dimension(:) :: csf
        integer i

        csf_spin = 0
        do i=1,size(csf)
            csf_spin = csf_spin - (real(csf(i))-1.5)
        enddo
    end function

    ! Calculates the number of possible S values (not their degeneracy)
    ! given a certain number of open shell electrons
    integer function num_S (nopen)
        integer, intent(in) :: nopen
        num_S = (nopen+2)/2
    end function

    ! Calculate the coefficients for each determinant contained in the
    ! CSF. These are calculated as the product of Clebsch-Gordon coeffs.
    ! working through the tree electron-by-electron. Each coeff. depends
    ! on the current total spin in the csf, the current total spin in
    ! the determinant and the spin of the current e-/posn being considered
    ! in either the determinant or the csf.
    ! dorder = the ordered list of alpha/beta for each spin orbital in det.
    real*8 pure function csf_coeff (csf, dorder, nopen)
        integer, intent(in), dimension(:) :: csf, dorder
        integer, intent(in) :: nopen
        real*8 S, M, scur, mcur, clb
        integer i

        S=0
        M=0
        csf_coeff = 1
        do i=1,nopen
            scur = -(csf(i)-1.5)
            mcur = -(real(dorder(i))-0.5)
            S = S + scur
            M = M + mcur

            clb = clbgrdn(S, M, scur, mcur)
            csf_coeff = csf_coeff * clb
            if (clb == 0) exit            
        enddo
    end function

    ! The Clebsch-Gordon coefficient with total spin S, projected
    ! component of spin M, with state with spin change 
    ! spin=(+/-)0.5 and projected spin change sig=(+/-)0.5
    real*8 pure function clbgrdn(S,M,spin,sig)
        real*8, intent(in) :: S,M,spin,sig
        !if (abs(spin).ne.0.5) then
        !    call stop_all ("clbgrdn","Spin incorrect")
        !endif
        if (spin.gt.0) then
            clbgrdn = sqrt((S+2*sig*M)/(2*S))
        else
            clbgrdn = -2*sig*sqrt((S+1-2*sig*M)/(2*(S+1)))
        endif
    end function

    ! Write a Yamanouchi symbol to output specified by nunit.
    ! if lTerm==.true. then terminate the line at the end.
    subroutine write_yama (nunit, yama, lTerm)
        integer, intent(in) :: nunit
        integer, intent(in), dimension(:) :: yama
        logical, intent(in) :: lTerm
        integer i

        write (nunit,'("(")',advance='no')
        do i=1,size(yama)
            write(nunit,'(i1)',advance='no') yama(i)
        enddo
        write(nunit,'(")")',advance='no')
        if (lTerm) write(nunit,*)
    end subroutine
end module

