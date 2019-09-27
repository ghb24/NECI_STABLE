! Module to collect and average the CI coefficients
! from equilibration until the end of the calculation
module sdt_amplitudes

  use bit_reps, only: extract_sign, decode_bit_det, encode_sign, niftot, nifd
  use constants, only: dp, lenof_sign, EPS, n_int, bits_n_int
  use DetBitOps, only: get_bit_excitmat, FindBitExcitLevel, EncodeBitDet
  use FciMCData, only: TotWalkers, iLutRef, CurrentDets, AllNoatHf, projedet, &
                       HashIndex, CurrentDets, ll_node, norm_psi, ilutHF, &
                       VaryShiftIter, iter
  use hash, only: hash_table_lookup,add_hash_table_entry,init_hash_table, &
                  clear_hash_table
  use LoggingData, only: n_store_ci_level,sorting_way
  use SystemData, only: nel,nbasis
  
  implicit none

  integer(n_int), allocatable :: ciCoeff_storage(:,:)
  integer :: hash_table_ciCoeff_size, first_free_entry = 0
  type(ll_node), pointer :: hash_table_ciCoeff(:)
  integer :: iout = 6
  integer :: nCyc, storSize
  real(dp), allocatable :: ciCoeff_storage_S(:,:)
  real(dp), allocatable :: ciCoeff_storage_D(:,:,:,:)
  real(dp), allocatable :: ciCoeff_storage_T(:,:,:,:,:,:)

contains


  subroutine init_ciCoeff
    nCyc = 0
    first_free_entry = 0
    hash_table_ciCoeff_size = 50000
    storSize = nbasis
    call init_hash_table(hash_table_ciCoeff)
    allocate(hash_table_ciCoeff(hash_table_ciCoeff_size))
    allocate(ciCoeff_storage(0:NIfTot,hash_table_ciCoeff_size))
  end subroutine init_ciCoeff


  subroutine fin_ciCoeff
    call clear_hash_table(hash_table_ciCoeff)
    deallocate(hash_table_ciCoeff)
    deallocate(ciCoeff_storage)
  end subroutine fin_ciCoeff

  !it prints not-averaged CI coeffs collected directly from last iter
  subroutine print_snapshot_ci_coeff
    integer :: i, ic, ex(2,3),icI
    real(dp) :: sign_tmp(lenof_sign)
    logical  :: tPar

    open (unit=21,file='SINGLES',status='replace')
    open (unit=22,file='DOUBLES',status='replace')
    open (unit=23,file='TRIPLES',status='replace')

    !main loop over the excitation level of the coeffs to be collected,
    !where 0 is the reference, 1 singles, 2 doubles and so on...
    do icI = 0, n_store_ci_level
      do i = 1, TotWalkers
        ic = 4
        call get_bit_excitmat(iLutRef(:,1),CurrentDets(:,i),ex,ic)
        call extract_sign(CurrentDets(:,i), sign_tmp)
        call GetBitExcitation(ilutRef(:,1),CurrentDets(:,i),ex,tPar)
        if(tPar) sign_tmp = -sign_tmp
        if(icI.eq.ic) then
          select case(icI)
          case(0)
            AllNoatHf = -sign_tmp
          case(1)
            write(21,'(G20.12,2I5)') sign_tmp/AllNoatHf(1),ex(1,1),ex(2,1)-nel
          case(2)
            write(22,'(G20.12,4I5)') sign_tmp/AllNoatHf(1),ex(1,1),ex(2,1)-nel,&
                                     ex(1,2),ex(2,2)-nel
          case(3)
            write(23,'(G20.12,6I5)') sign_tmp/AllNoatHf(1),ex(1,1),ex(2,1)-nel,&
                                     ex(1,2),ex(2,2)-nel,ex(1,3),ex(2,3)-nel
          end select
        end if
      end do
    end do

    close(21)
    close(22)
    close(23)
  end subroutine print_snapshot_ci_coeff


  !it prints averaged CI coeffs collected during the calcualtion
  subroutine print_averaged_ci_coeff
    integer :: i, ic, ex(2,3),icI,RefDet(nel)!!,nIEx(nel)
    real(dp) :: sign_tmp(lenof_sign)
    logical  :: tPar

    ! Kai told me there is already a subroutine to call in neci to sort
    ! the excitations, sth like: call sort(ciCoeff_storage), somewhere
    ! in quickLib?
    ! In the end I wrote my own sorting subroutine to organize the
    ! coefficients as needed depending on the case

    open (unit=31,file='SINGLES-AV',status='replace')
    open (unit=32,file='DOUBLES-AV',status='replace')
    open (unit=33,file='TRIPLES-AV',status='replace')

    write(iout,*) ''
!    write(iout,*) '-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --'
    write(iout,*) '-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --'
    write(iout,*) '*** CI COEFFICIENTS ***'
    write(iout,*) ''
    write(iout,*) 'Maximum excitation level of the CI coeffs =', n_store_ci_level

    do icI = 0, n_store_ci_level
       do i = 1, first_free_entry
          ic = 4
          call get_bit_excitmat(iLutRef(:,1),ciCoeff_storage(:,i),ex,ic)
!          ic = FindBitExcitLevel(ilutRef(:,1),ciCoeff_storage(:,i))
          if(icI.eq.ic) then
            call extract_sign(ciCoeff_storage(:,i),sign_tmp)
            call GetBitExcitation(ilutRef(:,1),ciCoeff_storage(:,i),ex,tPar)
            if(tPar) sign_tmp = -sign_tmp
             select case(icI)
             case(0)
                AllNoatHf = -sign_tmp/nCyc
                RefDet = projEDet(:,1)
!                write(iout,*) 'Total S+D+T = ', first_free_entry, 'nCyc=', nCyc
!                write(iout,*) sign_tmp/nCyc,AllNoatHf,sign_tmp(1),AllNoatHf(1),&
!                     (sign_tmp(1)/nCyc)/AllNoatHf(1)
             case(1)
                write(31,'(G20.12,2I5)') sign_tmp/(AllNoatHf(1)*nCyc), & 
                                         ex(1,1),ex(2,1)
             case(2) 
                write(32,'(G20.12,4I5)') sign_tmp/(AllNoatHf(1)*nCyc),ex(1,1),&
                                         ex(2,1),ex(1,2),ex(2,2)
             case(3)
                write(33,'(G20.12,6I5)') sign_tmp/(AllNoatHf(1)*nCyc),ex(1,1),&
                                         ex(2,1),ex(1,2),ex(2,2),ex(1,3),ex(2,3)
             end select
          end if
       end do
    end do

    close(31)
    close(32)
    close(33)
    call sorting(RefDet)

    write(iout,*) 'CI coefficients written in ASCII files'
    write(iout,*) ''
    write(iout,*) '-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --'

    call fin_ciCoeff()
  end subroutine print_averaged_ci_coeff

  
  subroutine storeCiCoeffs
    integer :: i, ic, ex(2,3),nIEx(nel)
    real(dp) :: sign_tmp(lenof_sign)
    nCyc = nCyc + 1

    ! loop through all occupied determinants 
    do i = 1, TotWalkers
       !definition of the maximal excitation level
       ic = 4
       !extraction of the excitation level from every determinant
       call get_bit_excitmat(iLutRef(:,1),CurrentDets(:,i),ex,ic)
       if(ic < 4) then
          call extract_sign(CurrentDets(:,i), sign_tmp)
          call decode_bit_det(nIEx, CurrentDets(:,i))
          call cache_sign(sign_tmp,nIEx,ex)
       end if
    end do

  end subroutine storeCiCoeffs


  !it updates the CI coeffs storage list
  subroutine cache_sign(sgn, nIEx, ex)
    integer :: hash_value ,ind
    integer, intent(in) :: nIEx(nel), ex(2,3)
    real(dp), intent(in) :: sgn(lenof_sign)
    real(dp) :: sign_tmp(lenof_sign)
    integer(n_int) :: ilut(0:NIfTot) 
    logical :: tSuccess

    ! encode the determinant into bit representation (ilut)
    call EncodeBitDet(nIEx,ilut)
    call hash_table_lookup(nIEx,ilut,NIfD,hash_table_ciCoeff,ciCoeff_storage,ind,hash_value,tSuccess)
    ! tSuccess is true when the coeff is found in the hash_table; so it gets updated
    if(tSuccess) then
       call extract_sign(ciCoeff_storage(:,ind),sign_tmp)
       sign_tmp = sign_tmp + sgn
!        if(all(nIEx.eq.projEDet(:,1))) then  !it is true only with the reference determinant
!           write(iout,*) 'sign_tmp = sign_tmp + sgn = ', sign_tmp, sgn
!        endif
    call encode_sign(ciCoeff_storage(:,ind),sign_tmp)

    ! tSuccess is false, then add a new entry to the CI coeffs storage
    else
       first_free_entry = first_free_entry + 1
       ! encode the determinant into bit representation (ilut)
       call EncodeBitDet(nIEx,ilut)
       ! store the encoded determinant in ciCoeff_storage
       ciCoeff_storage(:,first_free_entry) = ilut
       ! store the sign in ciCoeff_storage
       call encode_sign(ciCoeff_storage(:,first_free_entry),sgn)
       ! create a new hashtable entry
       call add_hash_table_entry(hash_table_ciCoeff,&
                                first_free_entry, hash_value)
    end if
  end subroutine cache_sign


  ! it allows to sort the averaged CI coeffs and list them in 2 different ways
  ! 1. OCC(alpha) + OCC(beta) + VIR(alpha) + VIR(beta) -> TODO
  ! 2. CLOSED orbs (alpha,beta) + OPEN orbs (alpha) + RELATIVE VIR orbs (beta) + LEFT VIR orbs (alpha,beta)
  subroutine sorting(RefDet)

    Implicit none
 
    double precision :: x
    double precision,allocatable :: S(:,:),D(:,:,:,:),T(:,:,:,:,:,:)
    integer :: i,j,k,a,b,c,z,p,Iopen(nel),iop,iMax,spo
    integer :: ial,ialMax,ialVir,ialVirMax,ibe,ibeMax,ibeVir,ibeVirMax    
    integer :: Ialpha(nel),Ibeta(nel),IalphaVir(nbasis-nel),IbetaVir(nbasis-nel)
    integer, intent(in) :: RefDet(nel)
    logical  :: check,noMatch,openEl,ClosedShellCase
    logical  :: aNoMatch,bNoMatch,cNoMatch,iNoMatch,jNoMatch,kNoMatch
 
    open (unit=101,file='SINGLES-AV',status='old', action='read')
    open (unit=102,file='DOUBLES-AV',status='old', action='read')
    open (unit=103,file='TRIPLES-AV',status='old', action='read')
    open (unit=31,file='SINGLES-AV-OR',status='replace')
    open (unit=32,file='DOUBLES-AV-OR',status='replace')
    open (unit=33,file='TRIPLES-AV-OR',status='replace')
 
    allocate(S(nel,nbasis))
    allocate(D(nel,nbasis,nel,nbasis))
    allocate(T(nel,nbasis,nel,nbasis,nel,nbasis))
    S(:,:)=0.d0
    D(:,:,:,:)=0.d0
    T(:,:,:,:,:,:)=0.d0

    openEl=.false.
    ClosedShellCase=.true.

    do z=2,nel
      p = RefDet(z) - RefDet(z-1)
      if((p.gt.1).and.(.not.openEl)) then
        Iopen(1)=RefDet(z-1)
        iop=2
        Iopen(iop)=RefDet(z)
        openEl=.true.
        ClosedShellCase=.false.
      else if((p.gt.1).and.(openEl)) then
        iop=iop+1
        Iopen(iop)=RefDet(z)
      endif
    enddo
    iMax=iop
    if (MOD(Iopen(iMax),2).eq.1) spo=1
    if (MOD(Iopen(iMax),2).eq.0) spo=-1


    ! normal reading for closed-shell systems
    if(ClosedShellCase) then

    write(iout,*) '***Reading coefficients in CLOSED-SHELL SYSTEM***'

      do
        read(101,*,IOSTAT=z) x,i,a
        if (z<0) then
          exit
        endif
        S(i,a)=x
      enddo
      close (101)

      do
        read(102,*,IOSTAT=z) x,i,a,j,b
        if (z<0) then
          exit
        endif
        D(i,a,j,b) = x
      enddo
      close (102)

      do
        read(103,*,IOSTAT=z) x,i,a,j,b,k,c
        if (z<0) then
          exit
        endif
        T(i,a,j,b,k,c) = x
      enddo
      close (103)


    ! open-shell systems require a specific reading in order to 
    ! reorganize the alpha and beta spin orbitals to make it
    ! readable in Molpro as efficiently as possible. Here there
    ! are implemented 2 different ways:
    else if(sorting_way.eq.1) then

    write(iout,*) 'Reading coefficients in OPEN-SHELL SYSTEM'
    write(iout,*) 'Coefficients listed in 1st way:'
    write(iout,*) '  OCC(alpha),OCC(beta),VIR(alpha),VIR(beta)'


      ial=0
      ibe=0
      do z=1,nel
        if (MOD(RefDet(z),2).eq.1) then
          ial=ial+1
          Ialpha(ial)=RefDet(z)
        write(iout,*) 'a_occ' ,ial,Ialpha(ial)
        else if (MOD(RefDet(z),2).eq.0) then
          ibe=ibe+1
          Ibeta(ibe)=RefDet(z)
        write(iout,*) 'b_occ' ,ibe,Ibeta(ibe)
        endif
      enddo
      ialMax=ial
      ibeMax=ibe

      ialVir=0
      ibeVir=0
      do i=1,nbasis
        j=0
        do z=1,nel
          if(i.eq.RefDet(z)) j=2
        enddo
        if(j.eq.2) cycle
        if (MOD(i,2).eq.1) then
          ialVir=ialVir+1
          IalphaVir(ialVir)=i
        write(iout,*) 'a' ,ialVir,IalphaVir(ialVir)
        else if (MOD(i,2).eq.0) then
          ibeVir=ibeVir+1
          IbetaVir(ibeVir)=i
        write(iout,*) 'b' ,ibeVir,IbetaVir(ibeVir)
        endif
      enddo
      ialVirMax=ialVir
      ibeVirMax=ibeVir

      write(iout,*) ialMax,ibeMax,ialVirMax,ibeVirMax


      do
        read(101,*,IOSTAT=z) x,i,a
        iNoMatch=.true.
        aNoMatch=.true.
        do ial=1,ialMax
          if(i.eq.Ialpha(ial).and.iNoMatch) then
            i=ial
            iNoMatch=.false.
          endif
        enddo
        do ibe=1,ibeMax
          if(i.eq.Ibeta(ibe).and.iNoMatch) then
            i=ialMax+ibe
            iNoMatch=.false.
          endif
        enddo
        do ialVir=1,ialVirMax
          if(a.eq.IalphaVir(ialVir).and.aNoMatch) then
            a=nel+ialVir
            aNoMatch=.false.
          endif
        enddo
        do ibeVir=1,ibeVirMax
          if(a.eq.IbetaVir(ibeVir).and.aNoMatch) then
            a=nel+ialVirMax+ibeVir
            aNoMatch=.false.
          endif
        enddo
        if (z<0) then
          exit
        endif
        S(i,a) = x
      enddo
      close (101)

      do
        read(102,*,IOSTAT=z) x,i,a,j,b
        iNoMatch=.true.
        jNoMatch=.true.
        aNoMatch=.true.
        bNoMatch=.true.
        do ial=1,ialMax
          if(i.eq.Ialpha(ial).and.iNoMatch) then
            i=ial
            iNoMatch=.false.
          else if(j.eq.Ialpha(ial).and.jNoMatch) then
            j=ial
            jNoMatch=.false.
          endif
        enddo
        do ibe=1,ibeMax
          if(i.eq.Ibeta(ibe).and.iNoMatch) then
            i=ialMax+ibe
            iNoMatch=.false.
          else if(j.eq.Ibeta(ibe).and.jNoMatch) then
            j=ialMax+ibe
            jNoMatch=.false.
          endif
        enddo
        do ialVir=1,ialVirMax
          if(a.eq.IalphaVir(ialVir).and.aNoMatch) then
            a=nel+ialVir
            aNoMatch=.false.
          else if(b.eq.IalphaVir(ialVir).and.bNoMatch) then
            b=nel+ialVir
            bNoMatch=.false.
          endif
        enddo
        do ibeVir=1,ibeVirMax
          if(a.eq.IbetaVir(ibeVir).and.aNoMatch) then
            a=nel+ialVirMax+ibeVir
            aNoMatch=.false.
          else if(b.eq.IbetaVir(ibeVir).and.bNoMatch) then
            b=nel+ialVirMax+ibeVir
            bNoMatch=.false.
          endif
        enddo
        if (z<0) then
          exit
        endif
        D(i,a,j,b) = x
      enddo
      close (102)

      do
        read(103,*,IOSTAT=z) x,i,a,j,b,k,c
        if (z<0) then
          exit
        endif
        T(i,a,j,b,k,c) = x
      enddo
      close (103)


    else if(sorting_way.eq.2) then
    ! TODO: in this way 2. of listing is not included yet the sorting of the open-shell
    !       RefDet with one or more virtual spatial orbs in between two occupied
    !       spin orbitals, e.g., p=4,6,8,... 

    write(iout,*) 'Reading coefficients in OPEN-SHELL SYSTEM'
    write(iout,*) 'Coefficients listed in 2nd way:'
    write(iout,*) '  CLOSED(alpha,beta),OPEN(alpha),OPEN(beta),VIRTUALS(alpha,beta)'

!      write(iout,*) 'PRINT OUTPUT: RefDet  = ', RefDet
!      write(iout,*) 'PRINT OUTPUT: p  = ', p
!      write(iout,*) 'PRINT OUTPUT: iMax = ', iMax
!      write(iout,*) 'PRINT OUTPUT: Iopen(1) = ', Iopen(1)
!      write(iout,*) 'PRINT OUTPUT: Iopen(2) = ', Iopen(2)
!      write(iout,*) 'PRINT OUTPUT: Iopen(3) = ', Iopen(3)
!      write(iout,*) 'PRINT OUTPUT: Iopen(iMax) = ', Iopen(iMax)
!      write(iout,*) 'PRINT OUTPUT: spo = ', spo
!      write(iout,*) 'PRINT OUTPUT: Iopen(iMax) = ', Iopen(iMax)
!      write(iout,*) 'PRINT OUTPUT: MOD(Iopen(iMax),2) = ', MOD(Iopen(iMax),2)
!      write(iout,*) 'PRINT OUTPUT: nel = ', nel
!      write(iout,*) 'PRINT OUTPUT: iop = ', iop
!      write(iout,*) 'PRINT OUTPUT: Iopen(iop)+spo = ', Iopen(iop)+spo
!      write(iout,*) 'PRINT OUTPUT: Iopen(1)+spo = ', Iopen(1)+spo
 

      do
        read(101,*,IOSTAT=z) x,i,a
        iNoMatch=.true.
        aNoMatch=.true.
        do iop=1,iMax
          if((i.eq.Iopen(iop)).and.(iNoMatch)) then
            i=nel-iMax+iop
            iNoMatch=.false.
          endif
          if((a.eq.Iopen(iop)+spo).and.(aNoMatch)) then
            a=nel+iop
            aNoMatch=.false.
          endif
        enddo
        S(i,a)=x
        if (z<0) then
          exit
        endif
      enddo
      close (101)
 
      do
        read(102,*,IOSTAT=z) x,i,a,j,b
        iNoMatch=.true.
        jNoMatch=.true.
        aNoMatch=.true.
        bNoMatch=.true.
        do iop=1,iMax
          if((i.eq.Iopen(iop)).and.(iNoMatch)) then
            i=nel-iMax+iop
            iNoMatch=.false.
          endif
          if((j.eq.Iopen(iop)).and.(jNoMatch)) then
            j=nel-iMax+iop
            jNoMatch=.false.
          endif
          if((a.eq.Iopen(iop)+spo).and.(aNoMatch)) then
!            write(iout,*) x,a, '-> nel+iop = ',nel, '+',iop,'=',nel+iop
            a=nel+iop
            aNoMatch=.false.
          endif
          if((b.eq.Iopen(iop)+spo).and.(bNoMatch)) then
            b=nel+iop
            bNoMatch=.false.
          endif
        enddo
          D(i,a,j,b) = x 
        if (z<0) then
          exit
        endif
      enddo
      close (102)
 
      do
        read(103,*,IOSTAT=z) x,i,a,j,b,k,c
        iNoMatch=.true.
        jNoMatch=.true.
        kNoMatch=.true.
        aNoMatch=.true.
        bNoMatch=.true.
        cNoMatch=.true.
        do iop=1,iMax
          if((i.eq.Iopen(iop)).and.(iNoMatch)) then
            i=nel-iMax+iop
            iNoMatch=.false.
          endif
          if((j.eq.Iopen(iop)).and.(jNoMatch)) then
            j=nel-iMax+iop
            jNoMatch=.false.
          endif
          if((k.eq.Iopen(iop)).and.(kNoMatch)) then
            k=nel-iMax+iop
            kNoMatch=.false.
          endif
          if((a.eq.Iopen(iop)+spo).and.(aNoMatch)) then
            a=nel+iop
            aNoMatch=.false.
          endif
          if((b.eq.Iopen(iop)+spo).and.(bNoMatch)) then
            b=nel+iop
            bNoMatch=.false.
          endif
          if((c.eq.Iopen(iop)+spo).and.(cNoMatch)) then
            c=nel+iop
            cNoMatch=.false.
          endif
        enddo
        T(i,a,j,b,k,c) = x
        if (z<0) then
          exit
        endif
      enddo
      close (103)

    endif


    do i = 1,nel
      do a = nel+1,nbasis 
        if(abs(S(i,a)).ne.0) then 
          write(31,'(G20.12,2I5)') S(i,a),i,a
!          write(31,'(G20.12,2I5)') S(i,a),i,a-nel
        end if
      end do
    end do

    do j = 2,nel
      do i = 1, nel-1
        do a = nel+1,nbasis 
          do b = a+1, nbasis
            if(abs(D(i,a,j,b)).ne.0 .and. i.lt.j) then 
              write(32,'(G20.12,4I5)') D(i,a,j,b),i,a,j,b
!              write(32,'(G20.12,4I5)') D(i,a,j,b),i,a-nel,j,b-nel
            end if
          end do
        end do
      end do
    end do

    do k = 3,nel
      do j = 2, nel-1
        do i = 1, nel-2
          do a = nel+1,nbasis 
            do b = a+1, nbasis
              do c = b+1, nbasis
                if(abs(T(i,a,j,b,k,c)).ne.0 .and. i.lt.j .and. j.lt.k) then 
                  write(33,'(G20.12,6I5)') T(i,a,j,b,k,c),i,a,j,b,k,c
!                  write(33,'(G20.12,6I5)') T(i,a,j,b,k,c),i,a-nel,j,b-nel,k,c-nel
                end if
              enddo
            enddo
          enddo
        enddo
      enddo
    end do
 
  end subroutine sorting

end module sdt_amplitudes
