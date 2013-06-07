! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
module gen_coul_ueg_mod
    use UMatCache, only: UMatInd, GTID
    use SystemData, only: BasisFN, nel, nBasisMax, nBasis, G1, NMSH, nMax, &
                          iSpinSkip, tHub, uHub, omega, iPeriodicDampingType,&
                          ALAT
    use IntegralsData, only: UMat, FCK
    use global_utilities
    use constants, only: dp, pi, pi2, THIRD
    implicit none

contains

    ! This is based on the gen_coul for electrons in a cubic box, but with
    ! a simpler coulomb integral
    subroutine gen_coul_hubnpbc ()

        ! This routine generates *ALL* possible combinations of Coulomb
        ! integrals stored in the form (u1 u2 | U | u1 u2) = UMAT(n1 n2 n3 n4)
        !
        ! This first call calculates the inner integral.
        ! The call to SCOUL calculates the outer integral.

        integer :: ll(3)
        integer :: a, b, ii, iSpinSkip
        real(dp) :: sums(3), sum
        integer :: i, j, k, l, id1, id2, id3, id4, ind
        type(timer), save :: proc_timer
        character(*), parameter :: this_routine = 'gen_coul_hubnpbc'

        proc_timer%timer_name = this_routine
        call set_timer(proc_timer)
        
        open (10, file='UMAT', status='unknown')
        ii = 0
        iSpinSkip = nBasisMax(2,3)
        ll(1) = nBasisMax(1,2) - nBasisMax(1,1) + 1
        ll(2) = nbasisMax(2,2) - nBasisMax(2,1) + 1
        ll(3) = nbasisMax(3,2) - nBasisMax(3,1) + 1
        do i = 1, nBasis, iSpinSkip
            do j = 1, nBasis, iSpinSkip
                do k = 1, nBasis, iSpinSkip
                    do l = 1, nBasis, iSpinSkip
                        sum = 1.0_dp
                        ! Original cal
                        ! call SLATCOULFOU (G1(1,I), G1(1,J), G1(1,K), &
                        !                   G1(1,L), NMSH, FCK, NMAX, SUM)
                        id1 = GTID(i)
                        id2 = GTID(j)
                        id3 = GTID(k)
                        id4 = GTID(l)
                        sum = uHub
                        do a = 1, 3
                            sums(a) = 0
                            do b = 1, ll(a)
                                sums(a) = sums(b) &
                                       + sin(b*pi*(G1(i)%K(a) - nBasisMax(a,1)+1) / (ll(a) + 1)) &
                                       * sin(b*pi*(G1(j)%K(a) - nBasismax(a,1)+1) / (ll(a) + 1)) &
                                       * sin(b*pi*(G1(k)%K(a) - nBasismax(a,1)+1) / (ll(a) + 1)) &
                                       * sin(b*pi*(G1(l)%K(a) - nBasismax(a,1)+1) / (ll(a) + 1))
                                       ! * sin(b*pi*G1(A,J)/(LL(A)+1))
                                       ! * sin(b*pi*G1(A,K)/(LL(A)+1))
                                       ! * sin(b*pi*G1(A,L)/(LL(A)+1))
                            enddo
                            sum = sum * sums(a) * (2.0_dp / (ll(a) + 1))**2
                        enddo

                        ! UMAT is stored as just spatial orbitals (not
                        ! spinorbitals). We store <IJ|KL>
                        ! Get the index of physical order UMAT element <IJ|KL>
                        ! Indices are internally reordered such that:
                        ! i >= k, j >= l, (i,k) >= (j,l)
                        UMAT(UMatInd(id1, id2, id3, id3, 0, 0)) = sum
                        if (abs(sum) > 1.0e-10_dp) &
                            write (10, '(4i7,f19.9)') id1, id2, id3, id4, sum
                    enddo
                enddo
            enddo
        enddo
        close (10)

        if (tHub) then
            ! V0 is subtracted from the diagonal element of the hamiltonian
            ! for the UEG, but we don't want this for the hubbard model.
            ! FCK is equiv to FCK(-nmsh/2:nmsh/2-1, ..., ...).
            ! Set FCK(0,0,0) = 0.
            ind = nbasis/2 + 1 + (nbasis/2) * nbasis * (nbasis + 1)
            FCK(ind) = (0, 0)
        endif

        call halt_timer (proc_timer)
    end subroutine


    ! THis is based upon the GEN_COUL for electrons in a cubic box but with a
    ! simpler coulomb integral
    subroutine gen_coul_ueg ()

        ! This routine generates *ALL* possible combinations of Coulomb 
        ! integrals stored in the form (u1 u2 | U | u1 u2) = UMAT(n1 n2 n3 n4)
        !
        ! This first call calculates the inner integral
        ! The call to SCOUL calculates the outer integral

        use sym_mod, only: roundsym, addelecsym, setupsym, lchksym

        type(BasisFn) :: ka, kb
        integer :: a, b, c, d, tx, ty
        integer :: ii, i, j, k, l, id1, id2, id3, id4, ind
        real(dp) :: lx, ly, p, q, t2, sum
        complex(dp) :: s, ci
        type(timer), save :: proc_timer
        character(*), parameter :: this_routine = 'gen_coul_ueg'

        proc_timer%timer_name = this_routine
        call set_timer (proc_timer)

        open (10, file='UMAT', status='unknown')
        ii = 0
        iSpinSkip = nBasismax(2,3)
        do i = 1, nBasis, iSpinSkip
            do j = 1, nBasis, iSpinSkip
                do k = 1, nBasis, iSpinSkip
                    do l = 1, nBasis, iSpinSkip
                        sum = 0.0_dp
                        ! Original cal
                        ! call SLATCOULFOU (G1(1,I), G1(1,J), G1(1,K), &
                        !                   G1(1,L), NMSH, FCK, NMAX, SUM)
                        id1 = GTID(i)
                        id2 = GTID(j)
                        id3 = GTID(k)
                        id4 = GTID(l)
                        if (tHub) then
                            call SetupSym (ka)
                            call SetupSym (kb)
                            call AddElecSym (k, G1, nBasisMax, ka)
                            call AddElecSym (l, G1, nBasisMax, ka)
                            call AddElecSym (i, G1, nBasisMax, kb)
                            call AddElecSym (j, G1, nBasisMax, kb)
                            call RoundSym (ka, nBasisMax)
                            call RoundSym (kb, nBasisMax)
                            !if ( (i == 1) .and. (j == 1) .and. (k <= l) ) &
                            !    write (40, '(6i3)', advance='no') &
                            !        (G1(ii,k), ii=1,2), (G1(ii,l), ii=1,2), &
                            !        (ka(ii), ii=1,2)
                            if (lChkSym(ka, kb)) then
                                !if (ChkMomEq (ka, kb, nBasisMax, 3)) then
                                    ! Omega for the hubbard model is the number 
                                    ! of sites in the lattice, (G1(ii,k),ii=1,2)
                                !    write (40,'(a10)',advance='no') "Y"
                                sum = uHub / omega
                            else
                                sum = 0
                            endif

                            if (.false.) then
                                s = 0.0_dp
                                ci = (0.0_dp, 1.0_dp)
                                ! Check the coulomb matrix elements by summing manually
                                c = G1(k)%K(1) + G1(l)%K(1) - G1(i)%K(1) - G1(j)%K(1)
                                d = G1(k)%K(2) + G1(l)%K(2) - G1(i)%K(2) - G1(j)%K(2)
                                tx = nBasisMax(1,4)
                                ty = nBasisMax(2,4)
                                t2 = (tx * tx) + (ty * ty)
                                lx = nBasisMax(1,5)
                                ly = nBasisMax(2,5)
                                do a = 1, nBasis, 2
                                    q = (G1(a)%K(1) * tx + G1(a)%K(2) * ty) / t2
                                    p = (G1(a)%K(1) * (-ty) + G1(a)%K(2) * tx) / t2
                                    s = s + exp(ci * (p*c/lx + q*d/ly) * 2 * pi) &
                                            * uHub / (omega**2)
                                enddo
                                write(41,"(2I3)",advance='no') c, d
                                if (abs(real(s) - sum) > 1.0e-7_dp) then
                                    write (41,*) sum, s
                                else
                                    write (41,*) sum, sum
                                endif
                                ! sum = real(s)
                            endif
                        else
                            ! We splice in the simple version of the fourier
                            ! transform
                            ! New, but not needed
                            ! call SetupSym (nBasisMax, ka)
                            ! call SetupSym (nBasisMax, kb)
                            ! call AddElecSym (k, G1, nBasisMax, ka)
                            ! call AddElecSym (l, G1, nBasisMax, ka)
                            ! call AddElecSym (i, G1, nBasisMax, kb)
                            ! call AddElecSym (j, G1, nBasisMax, kb)
                            ! call RoundSym (ka, nBasisMax)
                            ! call RoundSym (kb, nBasisMax)
                            !if (lChkSym(ka, kb)) then
                            a = G1(i)%K(1) - G1(k)%K(1)
                            b = G1(i)%K(2) - G1(k)%K(2)
                            c = G1(i)%K(3) - G1(k)%K(3)
                            if ( ((G1(l)%K(1) - G1(j)%K(1)) == a) .and. &
                                 ((G1(l)%K(2) - G1(j)%K(2)) == b) .and. &
                                 ((G1(l)%K(3) - G1(j)%K(3)) == c) .and. &
                                 ((a /= 0) .or. (b /= 0) .or. (c /= 0))) then
                                ! write (6,*) '(',i,j,'|',k,l,')',a,b,c
                                ! ajwt <ij|r_12^-1|kl> = v_(G_i-G_k) delta_((G_i-G_k)-(G_l-G_k))
                                ! v_G = 4 pi / G**2. G = 2 pi / L(nx, ny, nx) etc.
                                sum = ((a / ALAT(1))**2 + (b/ALAT(2))**2)
                                if (ALAT(3) /= 0.0_dp) sum = sum + (c / ALAT(3))**2
                                sum = 1 / (pi * sum * omega)
                            else
                                sum = 0.0_dp
                            endif

                        endif
                        ! UMAT is stored as just spatial orbitals (not
                        ! spin orbitals)
                        ! We store <IJ|KL>
                        ! Get the index of physical order UMAT element <IJ|KL>
                        ! Indices are internally ! reordered such that:
                        !i >= k, j >= l, (i,k) >= (j,l)
                        UMAT(UMatInd(id1, id2, id3, id4, 0, 0)) = sum
                        if (abs(sum) > 1.0e-10_dp) &
                            write (10, '(4i7,f19.9)') id1, id2, id3, id4, sum
                    enddo
                enddo
            enddo
        enddo
        close(10)
        if (tHub) then
            ! V0 is subtracted from the diagonal elements of the hamiltonian 
            ! for the UEG, but we don't want this for the hubbard model.
            ! FCK is equiv to FCK(-nmsh/2:nmsh/2-1, ..., ...).
            ! Set FCK(0,0,0) = 0.
            ind = nbasis/2 + 1 + (nbasis/2) * nbasis * (nbasis + 1)
            FCK(ind) = (0, 0)
        endif
        write (6, *)' !!! FINISHED CALCULATING ALL 2E INTEGRALS !!! '
  
        call halt_timer(proc_timer)

    end subroutine

    subroutine SlatCoulFouCou (G1, G2, G1P, G2P, N, CK, OUT)
    
        ! Returns the coulomb integral between the Slater determinants of
        ! plane wave basis, using the Fourier method.

        integer, intent(in) :: n
        integer, intent(in) :: G1(4), G2(4), G1P(4), G2P(4)
        complex(dp) :: CK(-N/2:N/2-1,-N/2:N/2-1,-N/2:N/2-1)
        real(dp), intent(out) :: OUT
        integer :: GD1, GD2, i, GD(3)
        real(dp) :: tot
        logical :: T, T2

        ! Check if the spins are correct
        if (G1(4) /= G1P(4) .or. G2(4) /= G2P(4)) then
            out = 0
            return
        endif

        T = .true.
        T2 = .true.
        tot = 0
        do i = 1, 3
            GD1 = G1(i) - G1P(i)
            GD2 = G2(i) - G2P(i)
            ! The 2* is to account for the fact that the UEG has basis fns
            ! which have k=2 pi n/L rather than pi n/L
            GD(i) = 2 * GD1
            if (GD1 /= GD2) T = .false.
            if (GD1 /= 0) T2 = .false.
        enddo
        
        if (T2) then
            ! GD1=0,0,0, which has to be removed as this cancels the
            ! background term.
            T = .false.
            ! write (6,*) CK(GD(1),GD(2),GD(3)),GD(1),GD(2),GD(3)
        endif
        ! if (T) then
            OUT = real(CK(2*GD(1),2*GD(2),2*GD(3)))
            ! OUT=4*3.1415926535_dp/(GD(1)*GD(1)+GD(2)*GD(2)+GD(3)*GD(3))
        !else
        !    OUT = 0.0_dp
        !endif
    end subroutine SlatCoulFouCou


    function ChkMomEq (k1, k2, nBasismax, kdim) result(lCmp)
        ! n.b. the third column of nBasisMax tells us whether it is tilted
        integer, intent(in) :: kDim, k1(kDim), k2(kDim), nBasisMax(5,*)
        integer :: j, lDim, kk1, kk2
        logical lCmp

        lCmp = .true.
        do j = 1, 3
            IF(NBASISMAX(1,3).EQ.1) THEN
               ! we have a tilted lattice
               IF(J.EQ.1) THEN
                  KK1=K1(1)+K1(2)
                  KK2=K2(1)+K2(2)
                  LDIM=NBASISMAX(1,2)*2
               ELSEIF(J.EQ.2) THEN
                  KK1=K1(1)-K1(2)
                  KK2=K2(1)-K2(2)
                  LDIM=NBASISMAX(1,2)*2
               ELSE
                  KK1=K1(J)
                  KK2=K2(J)
                  LDIM=NBASISMAX(J,2)-NBASISMAX(J,1)+1
               ENDIF

            ELSE
               KK1=K1(J)
               KK2=K2(J)
               LDIM=NBASISMAX(J,2)-NBASISMAX(J,1)+1
            ENDIF
            KK1=MOD(KK1,LDIM)
            IF(KK1.LT.NBASISMAX(J,1)) KK1=KK1+LDIM
            IF(KK1.GT.NBASISMAX(J,2)) KK1=KK1-LDIM
            KK2=MOD(KK2,LDIM)
            IF(KK2.LT.NBASISMAX(J,1)) KK2=KK2+LDIM
            IF(KK2.GT.NBASISMAX(J,2)) KK2=KK2-LDIM
            IF(KK1.NE.KK2) LCMP=.FALSE.
        enddo
        IF(KDIM.EQ.4.AND.K1(4).NE.K2(4)) LCMP=.FALSE.
    end function ChkMomEq
 
    function get_hub_umat_el (i, j, k, l) result(hel)
        use sym_mod, only: roundsym, addelecsym, setupsym, lchksym
        integer, intent(in) :: i, j, k, l
        HElement_t :: hel
        type(BasisFn) :: ka, kb

        call SetupSym (ka)
        call SetupSym (kb)
        call AddElecSym (k*2, G1, nBasisMax, ka)
        call AddElecSym (l*2, G1, nBasisMax, ka)
        call AddElecSym (i*2, G1, nBasisMax, kb)
        call AddElecSym (j*2, G1, nBasisMax, kb)
        call RoundSym (ka, nBasisMax)
        call RoundSym (kb, nBasisMax)
        if (lChkSym (ka, kb)) then
            hel = UMat (1)
        else
            hel = 0.0_dp
        endif
    end function

    function get_ueg_umat_el (idi, idj, idk, idl) result(hel)

        use SystemData, only: tUEG2, kvec, k_lattice_constant, dimen, Madelung
        integer, intent(in) :: idi, idj, idk, idl
        HElement_t :: hel
        integer :: i, j, k, l, a, b, c, iss, aneu
        real(dp) :: G, G2
        logical :: tCoulomb, tExchange          
        real(dp), parameter :: EulersConst = 0.5772156649015328606065120900824024_dp
        character(*), parameter :: this_routine = 'get_ueg_umat_el'

        !==================================================      
        if (tUEG2) then

            ! omit even numbers of idi if there is spin degeneracy
            ISS = nBasisMax(2,3) ! ick
            i = (idi - 1) * ISS + 1
            j = (idj - 1) * ISS + 1
            k = (idk - 1) * ISS + 1
            l = (idl - 1) * ISS + 1

            ! calucate unscaled momentum transfer
            a = kvec(i,1) - kvec(k, 1)
            b = kvec(i, 2) - kvec(k, 2)
            c = kvec(i, 3) - kvec(k, 3)

            ! Energy conservation
            if ((kvec(l, 1) - kvec(j, 1) == a) .and. &
            (kvec(l, 2) - kvec(j, 2) == b) .and. &
            (kvec(l, 3) - kvec(j, 3) == c) ) then

                ! no Coulomb (-> no divergency)
                if ( (a /= 0) .or. (b /= 0) .or. (c /= 0) ) then
                    ! Coulomb integrals are long-ranged, so calculated with 
                    ! 4 pi/G**2.
                    !AJWT <IJ|r_12^-1|KL> = v_(G_I-G_K) delta_((G_I-G_K)-(G_L-G_J)
                    ! v_G = 4 Pi/ G**2.  G=2 Pi/L(nx,ny,nx) etc.
                    ! For Coulomb interactions <ij|ij> we have explicitly excluded
                    ! the G=0 component as it is divergent.
                    ! This is the equivalent of adding a positive uniform 
                    ! background.
                    !scaled momentum transfer
                    G2 = (a *k_lattice_constant)**2 +(b *k_lattice_constant)**2 + (c *k_lattice_constant)**2         
                    ! check dimension
                    if(dimen == 3) then ! 3D
                        hel = (4.0_dp*PI) / (G2 * OMEGA)
                    else if (dimen ==2) then !2D
                        hel = (2.0_dp*PI) / (sqrt(G2) * OMEGA)
                    else if (dimen ==1) then !1D
                        hel = (-log(G2/4.0_dp) - 2.0_dp*EulersConst)/OMEGA
                    endif
                else  ! <ii|ii>
                    hel = 0
                endif  !Coulomb

            else  !no energy conservation
                hel = 0.0_dp
            endif
            return
        end if   !UEG2
!==================================================

        ISS = nBasisMax(2,3) ! ick

        i = (idi - 1) * ISS + 1
        j = (idj - 1) * ISS + 1
        k = (idk - 1) * ISS + 1
        l = (idl - 1) * ISS + 1

        tCoulomb = .false.
        tExchange = .false.

        if ( (i == k) .and. (j == l) ) tCoulomb = .true.
        if ( (i == l) .and. (j == k) ) tExchange = .true.

        ! The Uniform Electron Gas
        a = G1(i)%k(1) - G1(k)%k(1)
        b = G1(i)%k(2) - G1(k)%k(2)
        c = G1(i)%k(3) - G1(k)%k(3)

        if ( ((G1(l)%k(1) - G1(j)%k(1)) == a) .and. &
             ((G1(l)%k(2) - G1(j)%k(2)) == b) .and. &
             ((G1(l)%k(3) - G1(j)%k(3)) == c) ) then



            if ( (a /= 0) .or. (b /= 0) .or. (c /= 0) ) then
                ! WRITE(6,*) "(",I,J,"|",K,L,")",A,B,C
                ! Coulomb integrals are long-ranged, so calculated with 
                ! 4 pi/G**2.
                !AJWT <IJ|r_12^-1|KL> = v_(G_I-G_K) delta_((G_I-G_K)-(G_L-G_J)
                ! v_G = 4 Pi/ G**2.  G=2 Pi/L(nx,ny,nx) etc.
                ! For Coulomb interactions <ij|ij> we have explicitly excluded
                ! the G=0 component as it is divergent.
                ! This is the equivalent of adding a positive uniform 
                ! background.
                ! The effects
                G2 = ((a / ALAT(1))**2 +(b / ALAT(2))**2)
                if (ALAT(3) /= 0) G2 = G2 + (c / ALAT(3))**2

                ! Sum is G^2 / 4 Pi*Pi
                G = 2 * PI * sqrt(G2)
                hel = (1.0_dp / PI) / (G2 * ALAT(1) * ALAT(2) * ALAT(3))
 
                ! Sum is now (4 Pi / G^2 ) / Omega
                ! ALAT(4) is Rc, a cutoff length.
                ! For Exchange integrals we calculate <ij|kl>cell with a 
                ! potenial v(r)=1/r (r<Rc) and 0 (r>=Rc)
                if (tExchange) then
                    if (iPeriodicDampingType == 2) then
                        ! Spherical cutoff used.
                        ! For non-Coulomb integrals we calculate <ij|kl>_cell
                        ! with a potenial v(r)=1/r (r<Rc) and 0 (r>=Rc).
                        hel = hel * (1.0_dp - cos(G * ALAT(4)))
                    else if (iPeriodicDampingType == 1) then
                        ! Screened potential used.
                        ! For non-Coulomb integrals we calculate <ij|kl>
                        ! with a potenial v(r)=erfc(r/Rc)/r.
                        hel = hel * (1.0_dp - exp(-(G * ALAT(4))**2 / 2.0_dp))
                    endif
                endif
            else
                ! <ii|ii>
                hel = 0
                if (.not. tCoulomb .and. iPeriodicDampingType /= 0) then
                    ! The G=0 component is explicitly calculated for
                    ! non-Coulomb interactions as 2 PI Rc**2.
                    ! This is found by taking the Taylor expansion of the
                    ! attenuated and screened potentials and
                    ! considering the limit of G->0.  Both give the same
                    ! result.
                    call stop_all (this_routine, "<ii|ii> calculation error")
                    hel = 2 * pi * ALAT(4)**2 / (ALAT(1) * ALAT(2) * ALAT(3))
                endif
            endif
        else
            hel = 0
        endif

    end function

end module
