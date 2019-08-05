#include "macros.h" 

module gen_coul_ueg_mod
    use UMatCache, only: UMatInd, GTID
    use SystemData, only: BasisFN, nel, nBasisMax, nBasis, G1, NMSH, nMax, &
                          iSpinSkip, tHub, uHub, omega, iPeriodicDampingType,&
                          btHub, momIndexTable, breathingCont, tmodHub,&
                          ALAT,OrbECutoff,t_ueg_transcorr,t_ueg_3_body,&
                          tTranscorr,tContact,trpa_tc,tInfSumTCCalc,&
                          tInfSumTCPrint,tInfSumTCRead, Tperiodicinmom
    use IntegralsData, only: UMat, FCK
    use global_utilities
    use constants, only: sp, dp, pi, pi2, THIRD
    use Parallel_neci, only: iProcIndex,root
    use iso_c_hack
    use breathing_Hub, only: bHubIndexFunction
    implicit none
     
    interface
     function uu_tc_t (k2) result(u_tc)
        use constants, only: dp
        real(dp), intent(in) :: k2
        real(dp)  :: u_tc
     end function
    end interface
   
    procedure(uu_tc_t), pointer :: uu_tc

    save
! approximate the 3 body potential by 2 body terms and store them in  UMAT
     HElement_t(dp),  ALLOCATABLE :: UMAT_TC2(:,:,:)  
     HElement_t(dp),  ALLOCATABLE ::  UMAT_TC2_Contact(:), UMAT_TC2_Contact3D(:,:,:)
     real(dp) :: ktc_cutoff2,omega_p

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
                        UMAT(UMatInd(id1, id2, id3, id3)) = sum
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
            FCK(ind) = (0.0_sp, 0.0_sp)
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
                        UMAT(UMatInd(id1, id2, id3, id4)) = sum
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
            FCK(ind) = (0.0_sp, 0.0_sp)
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
        use constants, only: pi
        use SystemData, only: breathingCont
        HElement_t(dp) :: hel
        integer, intent(in) :: i, j, k, l
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
            if(tmodHub) hel = hel - breathingCont(bHubIndexFunction(i,j,k,l))
        else
            hel = 0.0_dp
        endif
    end function

      
!  the following fun is modified for transcorrelated Hamiltonian under RPA approx     
    function get_ueg_umat_el (idi, idj, idk, idl) result(hel)

        use SystemData, only: tUEG2, kvec, k_lattice_constant, dimen,PotentialStrength,TranscorrCutoff,nOccAlpha,nOccBeta
        integer, intent(in) :: idi, idj, idk, idl
        HElement_t(dp) :: hel
        integer :: i, j, k, l, a, b, c, iss, id1, id2, id3,id4,kmax,nsigma,sumind
        real(dp) :: G, G2,prefack,sprod
        real(dp) :: k_tc(3), pq_tc(3), u_tc, gamma_RPA, gamma_kmax,kveclength
        logical :: tCoulomb, tExchange          
        real(dp), parameter :: EulersConst = 0.5772156649015328606065120900824024_dp
        character(*), parameter :: this_routine = 'get_ueg_umat_el'


        ! Initialisation to satisfy compiler warnings
        hel = 0

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

        
     if(dimen==3)then
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
                

              if(t_ueg_transcorr)then  
! ============== the transcorrelated H under RPA =============================
                k_tc(1) = -2 * PI * a / ALAT(1)
                k_tc(2) = -2 * PI * b / ALAT(2)                
                k_tc(3) = -2 * PI * c / ALAT(3)
                G2 = k_tc (1) * k_tc (1) + k_tc (2) * k_tc (2) + k_tc (3) * k_tc (3)
                G = sqrt(G2)
                pq_tc(1) = (G1(k)%k(1) - G1(l)%k(1))*2*PI/ ALAT(1)
                pq_tc(2) = (G1(k)%k(2) - G1(l)%k(2))*2*PI/ ALAT(2)
                pq_tc(3) = (G1(k)%k(3) - G1(l)%k(3))*2*PI/ ALAT(3)
                
                 u_tc = uu_tc(G2)
                
                 if(t_ueg_3_body)then
                   hel = 4 * PI / G2 + G2 * u_tc - (pq_tc(1) * k_tc(1) +pq_tc(2) * k_tc(2) +pq_tc(3) * k_tc(3)) *u_tc  &
                         - UMAT_TC2(-a,-b,-c)
                 else
                   hel = 4 * PI / G2 + G2 * u_tc - (pq_tc(1) * k_tc(1) +pq_tc(2) * k_tc(2) +pq_tc(3) * k_tc(3)) *u_tc  &
                          - (nel-2) * G2 / (ALAT(1) * ALAT(2) * ALAT(3)) * u_tc**2&
                         - UMAT_TC2(-a,-b,-c)
                 end if       
               else
! =============== The original coulomb potential  ===================================              
                 G2 = ((a / ALAT(1))**2 +(b / ALAT(2))**2)
                 if (ALAT(3) /= 0) G2 = G2 + (c / ALAT(3))**2
 
                 ! Sum is G^2 / 4 Pi*Pi
                 G = 2 * PI * sqrt(G2)
                 hel = (1.0_dp / PI) / G2 
               end if 
                        
                 hel=hel / (ALAT(1) * ALAT(2) * ALAT(3)) 

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
                 hel = 0.
!                TC energy
                 if(t_ueg_transcorr)then
!                  if(i==j)then
!                   hel = -UMAT_TC2(0,0,0) / (ALAT(1) * ALAT(2) * ALAT(3))/4
!                    hel=0.0
!                  else 
!                   hel = -UMAT_TC2(0,0,0) / (ALAT(1) * ALAT(2) * ALAT(3))
!                  end if 
                  hel = -UMAT_TC2(0,0,0) / (ALAT(1) * ALAT(2) * ALAT(3))
                 end if 
                
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
     else if(dimen==2)then
        ! The Uniform Electron Gas
        a = G1(i)%k(1) - G1(k)%k(1)
        b = G1(i)%k(2) - G1(k)%k(2)

        if ( ((G1(l)%k(1) - G1(j)%k(1)) == a) .and. &
             ((G1(l)%k(2) - G1(j)%k(2)) == b) ) then



            if ( (a /= 0) .or. (b /= 0) ) then

              if(t_ueg_transcorr)then  
! ============== the transcorrelated H under RPA =============================
                k_tc(1) = -2 * PI * a / ALAT(1)
                k_tc(2) = -2 * PI * b / ALAT(2)                
                G2 = k_tc (1) * k_tc (1) + k_tc (2) * k_tc (2)
                G = sqrt(G2)
                pq_tc(1) = (G1(k)%k(1) - G1(l)%k(1))*2*PI/ ALAT(1)
                pq_tc(2) = (G1(k)%k(2) - G1(l)%k(2))*2*PI/ ALAT(2)
                
                if(G2<=ktc_cutoff2*(1.0+1.d-10))then
                 u_tc = 0.0
                else 
                 u_tc = - 2 * PI / G2/G
                end if
                
                
                 hel = 2 * PI / G + G2 * u_tc - (pq_tc(1) * k_tc(1) +pq_tc(2) * k_tc(2)) *u_tc  &
                        - (nel-2) * G2 / (ALAT(1) * ALAT(2)) * u_tc**2&
                        - UMAT_TC2(-a,-b,0)
               else
! =============== The original coulomb potential  ===================================              
                 G2 = ((a / ALAT(1))**2 +(b / ALAT(2))**2)
 
                 G = 2 * PI * sqrt(G2)
                 
                 hel = (2 * PI) / G 
               end if 
                        
                 hel=hel / (ALAT(1) * ALAT(2) ) 

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
                 hel = 0.
!                TC energy
                 if(t_ueg_transcorr)then
!                  if(i==j)then
!                    hel=0.
!                  else 
                   hel = -UMAT_TC2(0,0,0) / (ALAT(1) * ALAT(2))
!                  end if 
                 end if 
                
                if (.not. tCoulomb .and. iPeriodicDampingType /= 0) then
                    ! The G=0 component is explicitly calculated for
                    ! non-Coulomb interactions as 2 PI Rc**2.
                    ! This is found by taking the Taylor expansion of the
                    ! attenuated and screened potentials and
                    ! considering the limit of G->0.  Both give the same
                    ! result.
                    call stop_all (this_routine, "<ii|ii> calculation error")
                    hel = 2 * pi * ALAT(4)**2 / (ALAT(1) * ALAT(2))
                endif
            endif
        else
            hel = 0
        endif
      
      else
        write(6,*)'dimension error in get_ueg_umat_el',dimen
        stop
      end if  
      

    end function


!  the following fun is modified for three-particle excitations when one of the
!  creation and annhillaiton operator is the same from the minorirty component, and we sum-up all the
!  occupied orbitals
    function get_contact_umat_el_3b_sp(idi, idj, idk, idl) result(hel)

        use SystemData, only: kvec,dimen,PotentialStrength,TranscorrCutoff,nOccAlpha,nOccBeta,TranscorrGaussCutoff,t_trcorr_gausscutoff
        integer, intent(in) :: idi, idj, idk, idl
        HElement_t(dp) :: hel
        integer :: i, j, k, l, a, b, c, kmax,nsigma
        real(dp) ::  G2,prefack,GaussCutoff,GaussCutoff2,Gaussfact
        real(dp) :: kveclength
        character(*), parameter :: this_routine = 'get_contact_umat_el_3b'


        ! Initialisation to satisfy compiler warnings
        hel = 0

        !==================================================      
        i = idi
        j = idj
        k = idk
        l = idl

        if(dimen==3) then
        ! The Uniform Electron Gas
        a = G1(j)%k(1) - G1(l)%k(1)
        b = G1(j)%k(2) - G1(l)%k(2)
        c = G1(j)%k(3) - G1(l)%k(3)

        if ( ((G1(k)%k(1) - G1(i)%k(1)) == a) .and. &
             ((G1(k)%k(2) - G1(i)%k(2)) == b) .and. &
             ((G1(k)%k(3) - G1(i)%k(3)) == c) ) then

                        kmax=TranscorrCutoff
                        kveclength=dsqrt(dfloat(a**2+b**2+c**2))
                        prefack=2*PI/ALAT(1)
                        hel=0
                        if(kveclength.ge.kmax) then
                           G2 = kveclength*prefack
                           if(G1(i)%Ms.eq.-1) then
                                                nsigma=nOccAlpha
                           else
                                                nsigma=nOccBeta
                           endif
                         hel=-nsigma*4*PI**4/G2**4/ALAT(1)**3
                        endif

         endif !abc

        elseif(dimen==1) then

         a = G1(j)%k(1) - G1(l)%k(1)

         if ( G1(k)%k(1) - G1(i)%k(1) == a) then

            Gaussfact=1.d0
            if(t_trcorr_gausscutoff.and.a.ne.0) then
              GaussCutoff=TranscorrGaussCutoff*2.d0*PI/ALAT(1)
              GaussCutoff2=GaussCutoff**2
              Gaussfact=(1.d0-dexp(-GaussCutoff2*dfloat(a**2)))**2
            endif

                        kmax=TranscorrCutoff
                        kveclength=dabs(dfloat(a))
                        prefack=2*PI/ALAT(1)
                        hel=0
                        if(kveclength.ge.kmax.or.t_trcorr_gausscutoff) then
                           G2 = kveclength*prefack
                           if(G1(i)%Ms.eq.-1) then
                                                nsigma=nOccAlpha
                           else
                                                nsigma=nOccBeta
                           endif
                         hel=-nsigma*PotentialStrength**2/G2**2/ALAT(1)*Gaussfact
                        endif

          endif !a

        endif !dimen


   end function


!  the following fun is modified for three-particle excitations when one of the
!  creation and annhillaiton operator is the same from the minorirty component,
!  and we sum-up all the
!  occupied orbitals
    function get_contact_umat_el_3b_sap(idi, idj, idk, idl,nI) result(hel)

        use SystemData, only:kvec,dimen,PotentialStrength,TranscorrCutoff,nOccAlpha,nOccBeta,TranscorrGaussCutoff,t_trcorr_gausscutoff
        integer, intent(in) :: idi, idj, idk, idl, nI(nel)
        HElement_t(dp) :: hel
        integer :: i, j, k, l, a, b, c, sumind,sumind2,qocc,sumval
        integer :: alkvec(3), alexkvec(3), bekvec(3), beexkvec(3), kkvec(3), qocckvec(3)
        integer :: diffvec(3), diff2vec(3)
        real(dp) ::  G2,prefac, kveclength, sprod, difflength, diff2length
        real(dp) :: GaussCutoff,GaussCutoff2,Gaussfact
        character(*), parameter :: this_routine = 'get_contact_umat_el_3b'


        ! Initialisation to satisfy compiler warnings
        hel = 0

        !==================================================      

        if(is_beta(idi)) then 
                i = idj
                j = idi
                k = idl
                l = idk

        else
                i = idi
                j = idj
                k = idk
                l = idl

        endif

        if(dimen==3) then

          alkvec(1:3)=G1(i)%k(1:3)
          bekvec(1:3)=G1(j)%k(1:3)
          alexkvec(1:3)=G1(k)%k(1:3)
          beexkvec(1:3)=G1(l)%k(1:3)
                
          kkvec(1:3) = beexkvec(1:3) - bekvec(1:3)

          if ( ((alkvec(1) - alexkvec(1)) == kkvec(1)) .and. &
             ((alkvec(2) - alexkvec(2)) == kkvec(2)) .and. &
             ((alkvec(3) - alexkvec(3)) == kkvec(3)) ) then

                        prefac=1/4.d0/(ALAT(1) * ALAT(2) * ALAT(3))**2
                        do sumind=1, nel
                          qocc=nI(sumind)
                          qocckvec(1:3)=G1(qocc)%k(1:3)

                          if(is_beta(qocc)) then
                            if(qocc.eq.j) cycle

                            diffvec(1:3)=qocckvec(1:3)-bekvec(1:3)
                            sumval=0
                            do sumind2=1,3
                                sumval=sumval+diffvec(sumind2)**2
                            enddo
                            difflength=dsqrt(dfloat(sumval))
                            if(difflength.lt.TranscorrCutoff) cycle


                            diff2vec(1:3)=qocckvec(1:3)-beexkvec(1:3)
                            sumval=0
                            do sumind2=1,3
                                sumval=sumval+diff2vec(sumind2)**2
                            enddo
                            diff2length=dsqrt(dfloat(sumval))
                            if(diff2length.lt.TranscorrCutoff) cycle
                          else
                            if(qocc.eq.i) cycle
                            diffvec(1:3)=qocckvec(1:3)-alkvec(1:3)
                            sumval=0
                            do sumind2=1,3
                                sumval=sumval+diffvec(sumind2)**2
                            enddo
                            difflength=dsqrt(dfloat(sumval))
                            if(difflength.lt.TranscorrCutoff) cycle


                            diff2vec(1:3)=qocckvec(1:3)-alexkvec(1:3)
                            sumval=0
                            do sumind2=1,3
                                sumval=sumval+diff2vec(sumind2)**2
                            enddo
                            diff2length=dsqrt(dfloat(sumval))
                            if(diff2length.lt.TranscorrCutoff) cycle
                                
                          endif

                          sprod=0.d0
                          do sumind2=1,3
                             sprod=sprod+diffvec(sumind2)*diff2vec(sumind2)   
                          enddo
                          hel=hel-prefac*sprod/(difflength*diff2length)**3
                        enddo !sumind

           endif !abc

        elseif(dimen==1) then

         if(t_trcorr_gausscutoff.and.kkvec(1).ne.0) then
               GaussCutoff=TranscorrGaussCutoff*2.d0*PI/ALAT(1)
               GaussCutoff2=GaussCutoff**2
         endif

          alkvec(1)=G1(i)%k(1)
          bekvec(1)=G1(j)%k(1)
          alexkvec(1)=G1(k)%k(1)
          beexkvec(1)=G1(l)%k(1)

          kkvec(1) = beexkvec(1) - bekvec(1)

          if ( (alkvec(1) - alexkvec(1)) == kkvec(1)) then

                        prefac=Potentialstrength**2/ALAT(1)/(4.d0*PI**2)
                        do sumind=1, nel
                          qocc=nI(sumind)
                          qocckvec(1)=G1(qocc)%k(1)

                          if(is_beta(qocc)) then
                            if(qocc.eq.j) cycle

                            diffvec(1)=qocckvec(1)-bekvec(1)
                            difflength=dabs(dfloat(diffvec(1)))
                            if(.not.t_trcorr_gausscutoff.and.difflength.lt.TranscorrCutoff) cycle
                            


                            diff2vec(1)=qocckvec(1)-beexkvec(1)
                            diff2length=dabs(dfloat(diff2vec(1)))
                            if(.not.t_trcorr_gausscutoff.and.diff2length.lt.TranscorrCutoff) cycle

                          else
                            if(qocc.eq.i) cycle
                            diffvec(1)=qocckvec(1)-alkvec(1)
                            difflength=dabs(dfloat(diffvec(1)))
                            if(.not.t_trcorr_gausscutoff.and.difflength.lt.TranscorrCutoff) cycle


                            diff2vec(1)=qocckvec(1)-alexkvec(1)
                            diff2length=dabs(dfloat(diff2vec(1)))
                            if(.not.t_trcorr_gausscutoff.and.diff2length.lt.TranscorrCutoff) cycle

                          endif

                
                        if(diffvec(1).ne.0.and.diff2vec(1).ne.0) then
                          Gaussfact=1.d0
                          if(t_trcorr_gausscutoff) then
                                Gaussfact=(1.d0-dexp(-GaussCutoff2*dfloat(diffvec(1)**2)))*(1.d0-dexp(-GaussCutoff2*dfloat(diff2vec(1)**2)))
                          endif

                          hel=hel-prefac/(diffvec(1)*diff2vec(1))*Gaussfact
                        endif !diffvec(1).ne.0.and.diff2vec(1).ne.0

                        enddo !sumind

           endif !a

        endif!dimen


   end function

    
    
!  the following fun is modified for transcorrelated Hamiltonian under RPA approx     
    function get_contact_umat_el (idi, idj, idk, idl) result(hel)

        use SystemData, only: kvec, dimen,PotentialStrength,TranscorrCutoff,nOccAlpha,nOccBeta,TranscorrGaussCutoff,t_trcorr_gausscutoff,Tperiodicinmom
        integer, intent(in) :: idi, idj, idk, idl
        HElement_t(dp) :: hel
        integer :: i, j, k, l, a, b, c, kmax,nsigma,sumind
        real(dp) :: G2,prefack,sprod
        real(dp) :: k_tc(3), pq_tc(3), u_tc,kveclength
        real(dp) :: Gaussfact,GaussCutoff,GaussCutoff2
        logical :: tparallel, tmomconserv
        character(*), parameter :: this_routine = 'get_contact_umat_el'


        ! Initialisation to satisfy compiler warnings
        hel = 0

        !==================================================      
        i = idi
        j = idj 
        k = idk
        l = idl

        if(G1(i)%Ms.eq.G1(j)%Ms) then
                tparallel=.true.
                if(.not.trpa_tc) then
                        hel=0
                        return
                endif
        else
                tparallel=.false.
        endif 

     if(dimen==3)then
        ! The Uniform Electron Gas
        a = G1(j)%k(1) - G1(l)%k(1)
        b = G1(j)%k(2) - G1(l)%k(2)
        c = G1(j)%k(3) - G1(l)%k(3)

        if ( ((G1(k)%k(1) - G1(i)%k(1)) == a) .and. &
             ((G1(k)%k(2) - G1(i)%k(2)) == b) .and. &
             ((G1(k)%k(3) - G1(i)%k(3)) == c) ) then

                if (t_ueg_transcorr) then !transcorr

                        kmax=TranscorrCutoff
                        kveclength=dsqrt(dfloat(a**2+b**2+c**2))
                        prefack=2*PI/ALAT(1)
                        if(tparallel) then !parallel spin case

                                hel=0
                                if(kveclength.ge.kmax) then
                                        G2 = kveclength*prefack
                                        if(G1(i)%Ms.eq.-1) then
                                                nsigma=nOccAlpha
                                        else
                                                nsigma=nOccBeta
                                        endif
                                        hel=-nsigma*2*PI**2/G2**4/ALAT(1)**3

                                endif



                        else !anti-parallel spin case


                                hel = UMAT_TC2_Contact3D(abs(a),abs(b),abs(c))!/ PI**2
                                if(kveclength.ge.kmax ) then
                                         k_tc(1) = prefack * a
                                         k_tc(2) = prefack * b
                                         k_tc(3) = prefack * c
                                         G2 = kveclength*prefack
                                         pq_tc(1) = (G1(k)%k(1) -G1(l)%k(1))*prefack
                                         pq_tc(2) = (G1(k)%k(2) -G1(l)%k(2))*prefack
                                         pq_tc(3) = (G1(k)%k(3) -G1(l)%k(3))*prefack

                                         sprod=0.d0
                                         do sumind=1,3
                                                sprod=sprod+pq_tc(sumind)*k_tc(sumind)
                                         enddo

                                         hel = hel + 2*PI**2/G2 -sprod*2*PI**2/G2**3
                                end if

                        endif

                        hel = hel/ALAT(1)**3

!                       write(6,*)'in hel'
!                       write(6,*)hel, PotentialStrength, ALAT(1)

                else    !renormalization

                     hel=-4.d0*PI/(2.442749d0*dfloat(2*abs(NBASISMAX(1,2))+1))

                endif !transcorr or renormalization
          
        else
            hel = 0
        endif

      elseif(dimen==1) then

        a = G1(j)%k(1) - G1(l)%k(1)

        if(TranscorrCutoff.gt.0) then
                kmax=TranscorrCutoff
        else
                kmax=abs(NBASISMAX(1,2))
        endif

        if(Tperiodicinmom) then
                tmomconserv=mod(G1(k)%k(1) - G1(i)%k(1)-a,nbasis/2) == 0
                if(abs(a).gt.abs(NBASISMAX(1,2))) then
                        if(a.gt.0) then
                                a=a-nbasis/2
                        else
                                a=a+nbasis/2
                        endif
                endif
        else
                tmomconserv=((G1(k)%k(1) - G1(i)%k(1)) == a)
        endif

        if ( tmomconserv) then
         
          Gaussfact=1.d0
          if(t_trcorr_gausscutoff.and.a.ne.0) then
            GaussCutoff=TranscorrGaussCutoff*2.d0*PI/ALAT(1)
            GaussCutoff2=GaussCutoff**2
            Gaussfact=(1.d0-dexp(-GaussCutoff2*dfloat(a**2)))
          endif

         if(tparallel) then

           hel=0

           if(t_ueg_transcorr) then
              if(abs(a).ge.kmax.or.(t_trcorr_gausscutoff.and.a.ne.0)) then
                k_tc(1) = 2 * PI * a / ALAT(1)
                G2 = k_tc (1) * k_tc (1)
                if(G1(i)%Ms.eq.-1) then
                        nsigma=nOccAlpha
                else
                        nsigma=nOccBeta
                endif
                hel=hel-nsigma*PotentialStrength**2/G2/ALAT(1)*Gaussfact**2
              endif
           endif


         else


           hel = PotentialStrength

           if(t_ueg_transcorr)then
                  hel = hel+UMAT_TC2_Contact(abs(a)) !/ PI**2

                  k_tc(1) = 2 * PI * a / ALAT(1)
                  G2 = k_tc (1) * k_tc (1)
                  if( abs(a).ge.kmax.or.(t_trcorr_gausscutoff.and.a.ne.0) ) then
                      pq_tc(1) = (G1(k)%k(1) - G1(l)%k(1))*2*PI/ ALAT(1)
                      u_tc = -PotentialStrength/G2 !/PI

                      hel=hel-(PotentialStrength+pq_tc(1)*k_tc(1)*u_tc)*Gaussfact

                  end if
           end if

         endif !spin-parallel

           hel = hel/ALAT(1)

        else
            hel = 0
        endif ! ((G1(k)%k(1) - G1(i)%k(1)) == a)
 
      else
        write(6,*)'dimension error in get_contact_umat_el',dimen
        stop
      end if  
    
   end function 
 
    
    
      subroutine GEN_Umat_TC_Contact
        use SystemData, only: dimen,PotentialStrength,TranscorrCutoff,TranscorrIntCutoff
        use SystemData, only: tUnitary,t_trcorr_gausscutoff,TranscorrGaussCutoff
!       use Determinants, only: FDet
        use sym_mod, only: roundsym, addelecsym, setupsym, lchksym
        type(BasisFn) :: ka, kb
        integer :: i,j,k,shifti,shiftj,shiftk,diffi,diffj,diffk,sprodi,sprodij
        real(dp) :: sprod,length,difflength,GaussCutoff,GaussCutoff2
        integer :: AllocateStatus,kmax,kmaxcutoff
        integer :: maxj,mink,maxk,maxshiftj,maxshiftk,twokmax,signk
        integer :: shiftmaxj, shiftmaxk
        integer :: TrCutoffRead, TrIntCutoffRead, KmaxRead
        type(timer), save :: proc_timer
        logical :: tfile_exists

        real(dp) :: summasumma,prefactk,restofsum
        character(*), parameter :: this_routine = 'GEN_Umat_TC_Contact'
        character(LEN=100) :: dummyword


        proc_timer%timer_name = this_routine
        call set_timer (proc_timer)


      kmax=abs(NBASISMAX(1,2))
      twokmax=2*kmax
      if(t_trcorr_gausscutoff) then
        GaussCutoff=TranscorrGaussCutoff*2.d0*PI/ALAT(1)
        GaussCutoff2=GaussCutoff**2
        !It only calcualtes the value of the Gauss function 
        !if it is smaller then 1.e-12.
        TranscorrCutoff=ceiling(5.25652176976/GaussCutoff)+1
      endif
      if(TranscorrCutoff.gt.0) then
                kmaxcutoff=TranscorrCutoff
      else
                kmaxcutoff=kmax
      endif

      if(dimen==1) then

        ALLOCATE ( UMAT_TC2_Contact(0:2*kmax), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory for UMAT_TC2 ***"
        prefactk=-ALAT(1)*PotentialStrength**2/2.d0/PI**2

!       k=0
        summasumma=PI**2/6.d0
        do i=1, kmaxcutoff-1
           summasumma=summasumma-1.d0/dfloat(i)**2
        enddo !i

        UMAT_TC2_Contact(0)=summasumma


!       k=!0
        UMAT_TC2_Contact(1) = 1.d0/dfloat(kmaxcutoff)
        do i=2, 2*kmax
           UMAT_TC2_Contact(i)=UMAT_TC2_Contact(i-1)+1.d0/dfloat(kmaxcutoff+i-1)
        enddo !i

        do i=1, 2*kmax
           UMAT_TC2_Contact(i)=UMAT_TC2_Contact(i)/i
        enddo !i

        if(kmaxcutoff.le.kmax) then

          do i=2*kmaxcutoff, 2*kmax
           summasumma=0.d0
           do j=kmaxcutoff, i-kmaxcutoff
                summasumma=summasumma+1.d0/(dfloat(j)*dfloat(i-j))
           enddo

           UMAT_TC2_Contact(i)=UMAT_TC2_Contact(i)-summasumma/2.d0

          enddo !i
        endif

         if(t_trcorr_gausscutoff) then
            do i=0, 2*kmax
                 summasumma=0.d0
                 do j=-kmaxcutoff-i+1,kmaxcutoff+i-1
                    if(j.ne.0.and.j.ne.i) then
                       k=i-j
                       if(abs(j).lt.kmaxcutoff.or.abs(k).lt.kmaxcutoff) summasumma=summasumma+(1.d0-dexp(-GaussCutoff2*dfloat(j**2)))*(1.d0-dexp(-GaussCutoff2*dfloat(k**2)))/dfloat(j*k)
                    endif
                 enddo !j
                 UMAT_TC2_Contact(i)=UMAT_TC2_Contact(i)-summasumma/2.d0
            enddo  
         endif 

         UMAT_TC2_Contact(0:2*kmax)=prefactk*UMAT_TC2_Contact(0:2*kmax)

      elseif(dimen.eq.3.and.tUnitary) then
        ALLOCATE ( UMAT_TC2_Contact3D(0:twokmax,0:twokmax,0:twokmax), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory for UMAT_TC2 ***"

        if(tInfSumTCPrint) then
                open(32,file='TranscorrInfSum',status='unknown')
                write(32,*)"TrcorrCutoff=", kmaxcutoff
                write(32,*)"TrcorrIntCutoff=", TranscorrIntCutoff
                write(32,*)"Kmax=", kmax
        endif


        if(tInfSumTCCalc) then
        prefactk=(ALAT(1)/4.d0)
        !using the integral approximation
        restofsum=-4*PI/TranscorrIntCutoff
!       As the value is symmetric by exchanging i,j,k directions. We consider
!       only the i>=j>=k cases.
        do shifti=0,twokmax
!          shiftmaxj=min(int(dsqrt(dfloat(twokmax**2-shifti**2))),shifti)
           shiftmaxj=shifti
        do shiftj=0,shiftmaxj
!          shiftmaxk=min(int(dsqrt(dfloat(twokmax**2-shifti**2-shiftj**2))),shiftj)
           shiftmaxk=shiftj
        do shiftk=0,shiftmaxk
           summasumma=restofsum
                do i=-TranscorrIntCutoff+1,TranscorrIntCutoff-1
                        diffi=shifti-i
                        sprodi=i*diffi
                        maxj=int(dsqrt(dfloat(TranscorrIntCutoff**2-i**2)))
                do j=-maxj,maxj
                        diffj=shiftj-j
                        sprodij=j*diffj+sprodi
                        maxk=int(dsqrt(dfloat(TranscorrIntCutoff**2-i**2-j**2)))
                do k=-maxk,maxk
                        length=dsqrt(dfloat(i**2+j**2+k**2))
                        if(length.ge.TranscorrIntCutoff.or.length.lt.kmaxcutoff) cycle
                        diffk=shiftk-k
                        difflength=dsqrt(dfloat(diffi**2+diffj**2+diffk**2))
!                        write(6,*)'diff', diffi,diffj,diffk
                        if(difflength.lt.kmaxcutoff) cycle
                         sprod=dfloat(k*diffk+sprodij)
                        summasumma=summasumma+sprod/(length**3*difflength**3)
!                       write(6,*)'summasumma',i,j,k,summasumma,sprod/(length**3*difflength**3)
                enddo !k
                enddo !j
                enddo !i 
!                       write(6,*)'summasumma',shifti,shiftj,shiftk,summasumma
!                      if(shiftk.eq.1) stop
          summasumma=prefactk*summasumma
          UMAT_TC2_Contact3D(shifti,shiftj,shiftk)=summasumma
          if(tInfSumTCPrint) write(32,*)shifti,shiftj,shiftk,summasumma
          if(shifti.gt.shiftj)  then
                UMAT_TC2_Contact3D(shiftj,shifti,shiftk)=summasumma
                UMAT_TC2_Contact3D(shiftk,shiftj,shifti)=summasumma
                if(shiftj.gt.shiftk) then
                        UMAT_TC2_Contact3D(shiftj,shiftk,shifti)=summasumma
                        UMAT_TC2_Contact3D(shiftk,shifti,shiftj)=summasumma
                        UMAT_TC2_Contact3D(shifti,shiftk,shiftj)=summasumma
                endif
          else
              if(shiftj.gt.shiftk) then
                UMAT_TC2_Contact3D(shifti,shiftk,shiftj)=summasumma
                UMAT_TC2_Contact3D(shiftk,shifti,shiftj)=summasumma
              endif
          endif
        enddo !shiftk
        enddo !shiftj
        enddo !shifti

        elseif(tInfSumTCRead) then
                INQUIRE(FILE="TranscorrInfSum", EXIST=tfile_exists)
                if(.not.tfile_exists) stop "TranscorrInfSum is cannot be found!"
                open(32,file='TranscorrInfSum',status='unknown')

                read(32,*)dummyword, trcutoffread
                if(trcutoffread.ne.kmaxcutoff) then
                        write(6,*) "From TranscorrInfSum: ", trcutoffread, " From the input file:", kmaxcutoff
                        stop  "The value of TrCutoff is different in TranscorrInfSum."
                endif

                read(32,*)dummyword, trintcutoffread
                if(trintcutoffread.ne.TranscorrIntCutoff) then
                        stop "The value of TranscorrIntCutoff is different in TranscorrInfSum."
                endif

                read(32,*)dummyword, kmaxread
                if(kmaxread.lt.kmax) then
                        stop "The value of kmax is larger than kmaxread in TranscorrInfSum."
                endif

                do
                read(32,*,end=32) shifti,shiftj,shiftk,summasumma
                if(shifti.le.twokmax.and.shiftj.le.twokmax.and.shiftk.le.twokmax) then
                        UMAT_TC2_Contact3D(shifti,shiftj,shiftk)=summasumma
                        if(shifti.gt.shiftj)  then
                                UMAT_TC2_Contact3D(shiftj,shifti,shiftk)=summasumma
                                UMAT_TC2_Contact3D(shiftk,shiftj,shifti)=summasumma
                                if(shiftj.gt.shiftk) then
                                        UMAT_TC2_Contact3D(shiftj,shiftk,shifti)=summasumma
                                        UMAT_TC2_Contact3D(shiftk,shifti,shiftj)=summasumma
                                        UMAT_TC2_Contact3D(shifti,shiftk,shiftj)=summasumma
                                endif
                        else
                                if(shiftj.gt.shiftk) then
                                        UMAT_TC2_Contact3D(shifti,shiftk,shiftj)=summasumma
                                        UMAT_TC2_Contact3D(shiftk,shifti,shiftj)=summasumma
                                endif
                        endif
                endif
                enddo
32              write(6,*)"The values of inifinite sums are readed from TranscorrInfSum."
        endif


        if(tInfSumTCPrint.or.tInfSumTCRead) close(32)
      else
!       write(6,*) 'dimension error in GEN_Umat_TC', dimen
!       stop
      endif

!     write(6,*)'UMAT_TC2_Contact3D'

!     do i=0, twokmax
!     do j=0, twokmax
!     do k=0, twokmax
!       write(6,*)i,j,k,UMAT_TC2_Contact3D(i,j,k)
!     enddo
!     enddo
!     enddo
!     stop

      call halt_timer(proc_timer)
    end subroutine
    
    
     subroutine GEN_Umat_TC
        use SystemData, only: dimen
!        use Determinants, only: FDet
        use sym_mod, only: roundsym, addelecsym, setupsym, lchksym
        type(BasisFn) :: ka, kb
        integer :: a, b, c, d, tx, ty, kmax2, kmax2_cut
        integer :: ii, i, j, k, l, id1, id2, id3, id4, ind, AllocateStatus
        real(dp) :: lx, ly, p, q, t2, sum
        complex(dp) :: s, ci
        type(timer), save :: proc_timer

        real(dp) :: k_tc(3), pq_tc(3), gamma_RPA, gamma_kmax, G2, G,k1_tc(3),k2_tc(3)
        logical :: tCoulomb, tExchange          
        character(*), parameter :: this_routine = 'GEN_Umat_TC'
        

! ============== parameter for the correlation factor for TC method =============================
!        gamma_RPA=(4 * PI * nel/ (ALAT(1) * ALAT(2) * ALAT(3)))**0.25
!        gamma_kmax=sqrt(((NBASISMAX(1,2))*2*PI/ALAT(1))**2+ &
!                   ((NBASISMAX(2,2))*2*PI/ALAT(2))**2+((NBASISMAX(3,2))*2*PI/ALAT(3))**2)  
         ktc_cutoff2=OrbECutoff*(2*PI/ALAT(1))**2
         omega_p=dsqrt(4*pi*nel/ (ALAT(1) * ALAT(2) * ALAT(3)))
         
         
         
         if(t_ueg_3_body)then
          uu_tc => uu_tc_interpl
         else
          uu_tc => uu_tc_trunc
         end if 
        
        
        i=1       
        if(i==0)then
        
         call Madelungterm
         stop
        end if 
        
        
        
        open(10,file='3Bstatus',status='unknown')
        write(10,*)' With 3B RPA term -_-'
        close(10)

        proc_timer%timer_name = this_routine
        call set_timer (proc_timer)
      

!!!!!!!! generate UMAT_TC2, the Fourier transformation of (D u)^2
        kmax2=2*abs(NBASISMAX(1,2))
        
        
        
     if(dimen==3) then
     
        ALLOCATE ( UMAT_TC2(-kmax2:kmax2,-kmax2:kmax2,-kmax2:kmax2), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory for UMAT_TC2 ***"
        
        
        kmax2_cut=100
        
        do i = -kmax2, kmax2
           do j = -kmax2, kmax2
               do k = -kmax2, kmax2
                 UMAT_TC2(i,j,k) = 0.0
                 do id1= -kmax2_cut, kmax2_cut
                  k1_tc(1)= 2 * PI * id1 / ALAT(1)
                  k2_tc(1)= 2 * PI * (i-id1) / ALAT(1)
                 do id2= -kmax2_cut, kmax2_cut
                  k1_tc(2)= 2 * PI * id2 / ALAT(2)
                  k2_tc(2)= 2 * PI * (j-id2) / ALAT(2)
                 do id3= -kmax2_cut, kmax2_cut
                  k1_tc(3)= 2 * PI * id3 / ALAT(3)
                  k2_tc(3)= 2 * PI * (k-id3) / ALAT(3)
                  UMAT_TC2(i,j,k) = UMAT_TC2(i,j,k) - uu_tc_prod (k1_tc,k2_tc)
                 end do
                 end do
                 end do
                 UMAT_TC2(i,j,k) = UMAT_TC2(i,j,k) / (ALAT(1) * ALAT(2) * ALAT(3)) 
              end do
           end do
        end do  
        
      else if(dimen==2)then
      
        ALLOCATE ( UMAT_TC2(-kmax2:kmax2,-kmax2:kmax2,0:0), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory for UMAT_TC2 ***"
        

      
        kmax2_cut=1000
      
         do i = -kmax2, kmax2
           do j = -kmax2, kmax2
                 UMAT_TC2(i,j,0) = 0.0
                 do id1= -kmax2_cut, kmax2_cut
                  k1_tc(1)= 2 * PI * id1 / ALAT(1)
                  k2_tc(1)= 2 * PI * (i-id1) / ALAT(1)
                 do id2= -kmax2_cut, kmax2_cut
                  k1_tc(2)= 2 * PI * id2 / ALAT(2)
                  k2_tc(2)= 2 * PI * (j-id2) / ALAT(2)
                  k1_tc(3)= 0.0_dp
                  k2_tc(3)= 0.0_dp
                  
                   UMAT_TC2(i,j,0) = UMAT_TC2(i,j,0) - uu_tc_prod (k1_tc,k2_tc)
                 end do
                 end do
                 UMAT_TC2(i,j,0) = UMAT_TC2(i,j,0) / (ALAT(1) * ALAT(2) ) 
           end do
        end do  
     
      
      else 
        write(6,*) 'dimension error in GEN_Umat_TC', dimen
        stop
      end if  
         
!       open(10,file='UMAT_TC2',status='unknown')
!        do i=0,2
!        do j=0,2
!        do k=0,2
!         write(10,*),i,j,k,UMAT_TC2(i,j,k)
!        end do
!        end do
!        end do
!       close(10)
         
      
      
        call halt_timer(proc_timer)
        
    end subroutine
    

    
    function uu_tc_interpl0 (k2) result(u_tc)

     !   use SystemData, only: tUEG2, kvec, k_lattice_constant, dimen
        real(dp), intent(in) :: k2
        real(dp) :: u_tc
        
        if( k2<1.d-12 )then
         u_tc=0.0
        else 
         u_tc=-4*PI/k2/(k2+omega_p)
        end if
        
    end function
    
    function uu_tc_interpl (k2) result(u_tc)

     !   use SystemData, only: tUEG2, kvec, k_lattice_constant, dimen
        real(dp), intent(in) :: k2
        real(dp) :: u_tc
        
        if( k2<1.d-12 )then
         u_tc=0.0
        else 
         u_tc=-4*PI/k2/(k2+omega_p*2)
        end if
        
    end function
    

    
    function uu_tc_trunc (k2) result(u_tc)

        use SystemData, only: dimen
        real(dp), intent(in) :: k2
        real(dp) :: u_tc
      if(dimen==3)then  
        if(k2<=ktc_cutoff2*(1.0+1.d-10))then
         u_tc = 0.0
        else 
         u_tc = - 12.566370614359173 /k2/k2
        end if 
      else if (dimen==2)then
        if(k2<=ktc_cutoff2*(1.0+1.d-10))then
         u_tc = 0.0
        else 
         u_tc = - 6.283185307179586 / k2/dsqrt(k2)
        end if 
      else
        write(6,*) 'dimension error in uu_tc_prod', dimen
        stop
      end if

        
    end function
    
! Here we calculate the product k1*k2*u(k1)*u(k2) for the transcorrelated method
    function uu_tc_prod (k1_tc,k2_tc) result(uu_prod)

        use SystemData, only: tUEG2, kvec, k_lattice_constant, dimen
        real(dp), intent(in) :: k1_tc(3), k2_tc(3)
        real(dp) :: uu_prod, k1,k2,u1,u2,k12
        
        
      if(dimen==3)then  
        k12=k1_tc(1)*k2_tc(1)+k1_tc(2)*k2_tc(2)+k1_tc(3)*k2_tc(3)
        k1=k1_tc(1)*k1_tc(1)+k1_tc(2)*k1_tc(2)+k1_tc(3)*k1_tc(3)
        k2=k2_tc(1)*k2_tc(1)+k2_tc(2)*k2_tc(2)+k2_tc(3)*k2_tc(3)
        u1=uu_tc(k1)
        u2=uu_tc(k2)
        
        uu_prod=k12*u1*u2
      else if (dimen==2)then
        k12=k1_tc(1)*k2_tc(1)+k1_tc(2)*k2_tc(2)
        k1=k1_tc(1)*k1_tc(1)+k1_tc(2)*k1_tc(2)
        k2=k2_tc(1)*k2_tc(1)+k2_tc(2)*k2_tc(2)
        u1=uu_tc(k1)
        u2=uu_tc(k2)
        uu_prod=k12*u1*u2
      
      else
        write(6,*) 'dimension error in uu_tc_prod', dimen
        stop
      end if
      
    end function    

    subroutine Madelungterm
    use SystemData, only: dimen
    
    real(dp) :: kapa, sum, a,r
    integer  :: i,j,k,ii,m_cut
    
    
    if(dimen==3)then  !======================
    
    kapa=sqrt(pi)/ALAT(1)
    m_cut=20
    
    sum=0.0
    
    do i=-m_cut,m_cut
    do j=-m_cut,m_cut
    do k=-m_cut,m_cut
     
        ii=i*i+j*j+k*k
a=dsqrt(ii*ALAT(1)*ALAT(1))
        if(ii>0)then
        sum=sum+exp(-pi*ii)/pi/ii*(ALAT(1))**2
        end if
        end do
        end do
        end do

        sum=(sum-pi/kapa/kapa)/(ALAT(1))**3


        do i=-m_cut,m_cut
        do j=-m_cut,m_cut
        do k=-m_cut,m_cut

        ii=i*i+j*j+k*k
        if(ii>0)then
r=sqrt(1.0*ii)*ALAT(1) 
        sum=sum+erfc(kapa*r)/r
        end if
        end do
        end do
        end do
sum=sum-2*kapa/sqrt(pi)

        else !============================ dimen=2

kapa=sqrt(pi)/ALAT(1)
        m_cut=20

        sum=0.0

        do i=-m_cut,m_cut
        do j=-m_cut,m_cut

        ii=i*i+j*j
        if(ii>0)then
a=sqrt(1.0*ii)*2*pi/ALAT(1)
        sum=sum+erfc(a/2/kapa)/a
        end if
        end do
        end do


        sum=(sum*2*pi-2*dsqrt(pi)/kapa)/(ALAT(1))**2


        do i=-m_cut,m_cut
        do j=-m_cut,m_cut

        ii=i*i+j*j
        if(ii>0)then
r=sqrt(1.0*ii)*ALAT(1) 
        sum=sum+erfc(kapa*r)/r
        end if
        end do
        end do
sum=sum-2*kapa/sqrt(pi)


        end if !========================= not yet for other dimensions



        open(10,file='MadelungTerm.dat', status='unknown')
        write(10,*)sum*Nel/2
close(10)

        end subroutine

        !   prepare FCIDUMP for UEG, Debug    
        subroutine prep_ueg_dump

        use SystemData, only: tUEG2, kvec, k_lattice_constant, dimen
        use util_mod, only: get_free_unit
        HElement_t(dp) :: hel
        integer :: i, j, k, l, a, b, c, iss, id1, id2, id3, id4,ms2,nelec,i0,norb,i_unit
integer :: l1,l2,l3,r1,r2,r3,k1(3),k2(3),k3(3)
        real(dp) :: G, G2,energy,ak(3),bk(3),ck(3),a2,b2,c2
        real(dp) :: k_tc(3), pq_tc(3), u_tc, gamma_RPA, gamma_kmax
        logical :: tCoulomb, tExchange  
        character(*), parameter :: this_routine = 'prep_ueg_dump'
        namelist /FCI/ NORB,nelec,MS2

        if(iProcIndex == root) then

        ms2=0
        nelec=nel 
        norb=nBasis/2
i_unit=get_free_unit() 
        open(i_unit,file='FCIDUMP',status='unknown')
write(i_unit,FCI)
        do id1=1,norb
        do id2=1,norb
        do id3=1,norb
        do id4=1,norb
hel=get_ueg_umat_el (id1, id2, id3, id4)
        if(abs(hel).gt.1.d-30)then
        write(i_unit,'(1x,e23.16,4I4)')hel,id1,id3,id2,id4
        end if
        end do
        end do
        end do
        end do
        i0=0
        do i=1,norb
        do j=i,i

        id1=(i-1)*2+1

        a = G1(id1)%k(1) 
        b = G1(id1)%k(2) 
c = G1(id1)%k(3)

        hel=a*a/2.0+b*b/2.0+c*c/2.0
        hel=hel*(2*pi/ALAT(1))**2

        if(abs(hel).gt.1.d-30)then
        write(i_unit,'(1x,e23.16,4I4)')hel,i,j,i0,i0
        end if
        end do
        end do
        energy=0.0
        write(i_unit,'(1x,e23.16,4I4)')energy,i0,i0,i0,i0
close(i_unit)

        !========= L mat ========================
        if(.not.t_ueg_3_body)stop
i_unit=get_free_unit() 
        open(i_unit,file='TCDUMP',status='unknown',access='append')
        do l1=1,norb
        do l2=l1,norb
        do l3=l2,norb
        do r1=1,norb
        do r2=2,norb
        do r3=3,norb

hel=get_lmat_ueg(l1,l2,l3,r1,r2,r3)
        if(dabs(hel)>1.d-8)then
        write(i_unit,'(1x,e23.16,6I4)')hel,l1,l2,l3,r1,r2,r3
        !           print '(1x,e23.16,6I4)', hel,l1,l2,l3,r1,r2,r3
        end if

        end do
        end do
        end do
        end do
        end do
        end do




close(i_unit)
        end if
        stop
        end subroutine 

function get_lmat_ueg (l1,l2,l3,r1,r2,r3) result(hel)

        use SystemData, only: dimen
        integer, intent(in) :: l1,l2,l3,r1,r2,r3
integer :: i,j,k,a,b,c,k1(3),k2(3),k3(3)
        real(dp) :: hel,ak(3),bk(3),ck(3),a2,b2,c2
        i=(l1-1)*2+1
        j=(l2-1)*2+1
        k=(l3-1)*2+1
        a=(r1-1)*2+1
        b=(r2-1)*2+1
        c=(r3-1)*2+1


        !============ copy from above

        if(dimen==3)then
        ! The Uniform Electron Gas
        k1(1) = G1(i)%k(1) - G1(a)%k(1)
        k1(2) = G1(i)%k(2) - G1(a)%k(2)
k1(3) = G1(i)%k(3) - G1(a)%k(3)

        k2(1) = G1(j)%k(1) - G1(b)%k(1)
        k2(2) = G1(j)%k(2) - G1(b)%k(2)
k2(3) = G1(j)%k(3) - G1(b)%k(3)

        k3(1) = G1(k)%k(1) - G1(c)%k(1)
        k3(2) = G1(k)%k(2) - G1(c)%k(2)
k3(3) = G1(k)%k(3) - G1(c)%k(3)

        if(all((k1+k2+k3)==0)) then
        ak(1) = -2 * PI * k1(1) / ALAT(1)
        ak(2) = -2 * PI * k1(2) / ALAT(2)                
        ak(3) = -2 * PI * k1(3) / ALAT(3)
        bk(1) = -2 * PI * k2(1) / ALAT(1)
        bk(2) = -2 * PI * k2(2) / ALAT(2)                
        bk(3) = -2 * PI * k2(3) / ALAT(3)
        ck(1) = -2 * PI * k3(1) / ALAT(1)
        ck(2) = -2 * PI * k3(2) / ALAT(2)                
ck(3) = -2 * PI * k3(3) / ALAT(3)


        a2=ak(1)*ak(1)+ak(2)*ak(2)+ak(3)*ak(3)
        b2=bk(1)*bk(1)+bk(2)*bk(2)+bk(3)*bk(3)
c2=ck(1)*ck(1)+ck(2)*ck(2)+ck(3)*ck(3)

        hel= uu_tc(a2)*uu_tc(b2)*(ak(1)*bk(1)+ak(2)*bk(2)+ak(3)*bk(3))&
        +uu_tc(a2)*uu_tc(c2)*(ak(1)*ck(1)+ak(2)*ck(2)+ak(3)*ck(3))&
+uu_tc(b2)*uu_tc(c2)*(bk(1)*ck(1)+bk(2)*ck(2)+bk(3)*ck(3))
        !           hel=hel/(ALAT(1) * ALAT(2) * ALAT(3))**2/3
        ! a permutation factor 3 is missing in other place, so here we ignore the /3
        hel=hel/(ALAT(1) * ALAT(2) * ALAT(3))**2

        else
        hel=0.0
        end if
        else
        print *, 'at moment Lmat is only available for 3D UEG'
        stop
        end if

        end function

function get_lmat_ua (l1,l2,l3,r1,r2,r3) result(hel)
        use SystemData, only: dimen, TranscorrCutoff,PotentialStrength,TranscorrGaussCutoff,t_trcorr_gausscutoff
        integer, intent(in) :: l1,l2,l3,r1,r2,r3
        integer :: i,j,k,a,b,c,k1(3),k2(3),k3(3)
        real(dp) :: hel,ak(3),bk(3),ck(3),a3,b3,c3,klength1, klength2
        real(dp) :: GaussCutoff,GaussCutoff2,Gaussfact
!       logical :: tSpinCorrect


     if(dimen==3)then

!       The spin has been set before in get_lmat_el_ua l1,l2,r1,r2 have the same
!       spin. l3 and r3 have the opposite spin.
        i=l1
        j=l2
        k=l3
        a=r1
        b=r2
        c=r3

        k1(1) = G1(i)%k(1) - G1(a)%k(1)
        k1(2) = G1(i)%k(2) - G1(a)%k(2)
        k1(3) = G1(i)%k(3) - G1(a)%k(3)

        klength1=dsqrt(dfloat(k1(1)**2+k1(2)**2+k1(3)**2))

        k2(1) = G1(b)%k(1) - G1(j)%k(1)
        k2(2) = G1(b)%k(2) - G1(j)%k(2)
        k2(3) = G1(b)%k(3) - G1(j)%k(3)

        klength2=dsqrt(dfloat(k2(1)**2+k2(2)**2+k2(3)**2))

        k3(1) = G1(c)%k(1) - G1(k)%k(1)
        k3(2) = G1(c)%k(2) - G1(k)%k(2)
        k3(3) = G1(c)%k(3) - G1(k)%k(3)

        if(all((k2-k1+k3)==0).and.klength1.ge.TranscorrCutoff.and.klength2.ge.TranscorrCutoff) then
           ak(1) = 2 * PI * k1(1) / ALAT(1)
           ak(2) = 2 * PI * k1(2) / ALAT(2)
           ak(3) = 2 * PI * k1(3) / ALAT(3)
           bk(1) = 2 * PI * k2(1) / ALAT(1)
           bk(2) = 2 * PI * k2(2) / ALAT(2)
           bk(3) = 2 * PI * k2(3) / ALAT(3)


           a3=dsqrt(ak(1)*ak(1)+ak(2)*ak(2)+ak(3)*ak(3))**(-3)
           b3=dsqrt(bk(1)*bk(1)+bk(2)*bk(2)+bk(3)*bk(3))**(-3)

           hel=-2.d0*PI**4*a3*b3*(ak(1)*bk(1)+ak(2)*bk(2)+ak(3)*bk(3)) !&

         else

               hel=0.d0
         end if

        hel=hel/(ALAT(1) * ALAT(2) * ALAT(3))**2
      elseif(dimen==1) then


        i=l1
        j=l2
        k=l3
        a=r1
        b=r2
        c=r3

        k1(1) = G1(i)%k(1) - G1(a)%k(1)

        klength1=dabs(dfloat(k1(1)))

        k2(1) = G1(b)%k(1) - G1(j)%k(1)

        klength2=dabs(dfloat(k2(1)))

        k3(1) = G1(c)%k(1) - G1(k)%k(1)

        if(k2(1)-k1(1)+k3(1)==0.and.((klength1.ge.TranscorrCutoff.and.klength2.ge.TranscorrCutoff).or.(t_trcorr_gausscutoff.and.(k1(1).ne.0.and.k2(1).ne.0)))) then

           Gaussfact=1.d0
           if(t_trcorr_gausscutoff) then
            GaussCutoff=TranscorrGaussCutoff*2.d0*PI/ALAT(1)
            GaussCutoff2=GaussCutoff**2
            Gaussfact=(1.d0-dexp(-GaussCutoff2*dfloat(k1(1)**2)))*(1.d0-dexp(-GaussCutoff2*dfloat(k2(1)**2)))
           endif

           ak(1) = 2 * PI * k1(1) / ALAT(1)
           bk(1) = 2 * PI * k2(1) / ALAT(1)


           hel=-PotentialStrength**2/(ak(1)*bk(1))/2.d0*Gaussfact

         else

               hel=0.d0
         end if
                
        hel=hel/ALAT(1)
      else
       print *, 'at moment Lmat is only available for 1D and 3D contact interactions'
       stop
      end if


    end function

end module
