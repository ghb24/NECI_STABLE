! Module to collect and average the CI coefficients
! from equilibration until the end of the calculation
module sdt_amplitudes

  use bit_reps, only: extract_sign, decode_bit_det, encode_sign, niftot, nifd
  use constants, only: dp, lenof_sign, EPS, n_int, bits_n_int,int64,iout
  use DetBitOps, only: get_bit_excitmat, FindBitExcitLevel, EncodeBitDet
  use util_mod, only: near_zero, operator(.isclose.)
  use FciMCData, only: TotWalkers, iLutRef, CurrentDets, AllNoatHf, projedet, &
                       HashIndex, CurrentDets, ll_node, norm_psi, ilutHF, &
                       VaryShiftIter, iter
  use hash, only: hash_table_lookup,add_hash_table_entry,init_hash_table, &
                  clear_hash_table
  use LoggingData, only: n_store_ci_level,sorting_way,n_iter_after_equ
  use SystemData, only: nel,nbasis
  use sort_mod, only: sort
  use DetBitOps, only: ilut_lt, ilut_gt

  implicit none

  integer(n_int), allocatable :: ciCoeff_storage(:,:)
  integer :: hash_table_ciCoeff_size, first_free_entry
  type(ll_node), pointer :: hash_table_ciCoeff(:)
  integer :: nCyc
  character (len=90) :: fileCICoeffSnpsht,fileCICoeffAv,fileCIcoeffSort
!  integer :: storSize
!  real(dp), allocatable :: ciCoeff_storage_S(:,:)
!  real(dp), allocatable :: ciCoeff_storage_D(:,:,:,:)
!  real(dp), allocatable :: ciCoeff_storage_T(:,:,:,:,:,:)

contains


  subroutine init_ciCoeff
    nCyc = 0
    first_free_entry = 0
    hash_table_ciCoeff_size = 50000
!    storSize = nbasis
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
    integer :: ic, ex(2,3),icI
    integer(int64) :: i
    real(dp) :: sign_tmp(lenof_sign)
    logical  :: tPar

    !main loop over the excitation level of the coeffs to be collected,
    !where 0 is the reference, 1 singles, 2 doubles and so on...
    do icI = 0, n_store_ci_level
       if(icI.ne.0) then
         write(fileCICoeffSnpsht, '( "CI_COEFF_",I1,"_snapshot" )') icI
         open (unit=20+icI,file=fileCICoeffSnpsht,status='replace')
       endif
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
      if(icI.ne.0) close(20+icI)
    end do

  end subroutine print_snapshot_ci_coeff


  !it prints averaged CI coeffs collected during the calcualtion
  subroutine print_averaged_ci_coeff
    integer :: i, ic, ex(2,3),icI,RefDet(nel)!!,nIEx(nel)
    real(dp) :: sign_tmp(lenof_sign)
    logical  :: tPar

    ! There is a subroutine to call in neci to sort the excitations,
    ! lib/quicksort.F90.template: the relevant subroutine is:
    ! sort_mod::sort
    ! In the end I wrote my own sorting subroutine to organize the
    ! coefficients as needed depending on the case

!    call sort(ciCoeff_storage(:,1:first_free_entry), ilut_lt, ilut_gt)
!    call sort(ciCoeff_storage(:,1:first_free_entry), indices_lt, indices_gt)

    write(iout,*) ''
!    write(iout,*) '-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --'
    write(iout,*) '-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --'
    write(iout,*) ''
    write(iout,*) '*** CI COEFFICIENTS ***'
      if (n_iter_after_equ.lt.0) then
        write(iout,*) ''
        write(iout,*) "!CI coefficients collection better starts after equilibration:"
        write(iout,*) "  -> set number of iterations after equilibration greater or equal to 0"
      endif
    write(iout,*) ''
    write(iout,*) 'Maximum excitation level of the CI coeffs =', n_store_ci_level

    do icI = 0, n_store_ci_level
       if(icI.ne.0) then
         write(fileCICoeffAv, '( "CI_COEFF_",I1,"_AV" )') icI
         open (unit=30+icI,file=fileCICoeffAv,status='replace')
       endif
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
       if(iCI.ne.0) close(30+icI)
    end do

    call sorting(RefDet)

    write(iout,*) '-> CI coefficients written in ASCII files'
    write(iout,*) ''
    write(iout,*) '-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --'

    call fin_ciCoeff()
  end subroutine print_averaged_ci_coeff


  subroutine storeCiCoeffs
    integer :: ic, ex(2,3),nIEx(nel)
    integer(int64) :: i
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
          call cache_sign(sign_tmp, nIEx)
       end if
    end do

  end subroutine storeCiCoeffs


  !it updates the CI coeffs storage list
  subroutine cache_sign(sgn, nIEx)
    integer, intent(in) :: nIEx(nel)
    real(dp), intent(in) :: sgn(lenof_sign)

    integer :: hash_value, ind
    real(dp) :: sign_tmp(lenof_sign)
    integer(n_int) :: ilut(0:NIfTot)
    logical :: tSuccess

    ! encode the determinant into bit representation (ilut)
    call EncodeBitDet(nIEx, ilut)
    call hash_table_lookup(nIEx, ilut, NIfD, hash_table_ciCoeff, &
                           ciCoeff_storage, ind, hash_value, tSuccess)
    ! tSuccess is true when the coeff is found in the hash_table; so it gets updated
    if(tSuccess) then
       call extract_sign(ciCoeff_storage(:, ind), sign_tmp)
       sign_tmp = sign_tmp + sgn
!        if(all(nIEx.eq.projEDet(:,1))) then  !it is true only with the reference determinant
!           write(iout,*) 'sign_tmp = sign_tmp + sgn = ', sign_tmp, sgn
!        endif
    call encode_sign(ciCoeff_storage(:, ind), sign_tmp)

    ! tSuccess is false, then add a new entry to the CI coeffs storage
    else
       first_free_entry = first_free_entry + 1
       ! encode the determinant into bit representation (ilut)
       call EncodeBitDet(nIEx, ilut)
       ! store the encoded determinant in ciCoeff_storage
       ciCoeff_storage(:, first_free_entry) = ilut
       ! store the sign in ciCoeff_storage
       call encode_sign(ciCoeff_storage(:, first_free_entry), sgn)
       ! create a new hashtable entry
       call add_hash_table_entry(hash_table_ciCoeff, &
                                 first_free_entry, hash_value)
    end if
  end subroutine cache_sign


  ! it sorts the averaged CI coeffs and list them in different ways.
  ! the sorting way can be chosen from the input for open-shell systems
  subroutine sorting(RefDet)

    Implicit none

    integer, intent(in) :: RefDet(nel)
    double precision :: x
    double precision,allocatable :: S(:,:),D(:,:,:,:),T(:,:,:,:,:,:)
    integer :: i,j,k,a,b,c,z,p,Iopen(nel),iop,iMax,spo,Nind
    integer :: ial,ialMax,ialVir,ialVirMax,ibe,ibeMax,ibeVir,ibeVirMax,Itot(nbasis),signCI
    integer ::  indCoeff(n_store_ci_level),indCoeffV(n_store_ci_level),indCoef(2,n_store_ci_level)
    logical  :: check,noMatch,openEl,ClosedShellCase
    logical  :: aNoMatch,bNoMatch,cNoMatch,iNoMatch,jNoMatch,kNoMatch
!    integer :: Ialpha(nel),Ibeta(nel),IalphaVir(nbasis-nel),IbetaVir(nbasis-nel)
!    double precision :: x,S(nel,nbasis-nel),D(nel,nbasis-nel,nel,nbasis-nel)
!    double precision :: T(nel,nbasis-nel,nel,nbasis-nel,nel,nbasis-nel)

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
      else if(p.eq.1.and.(.not.openEl).and.z.eq.nel.and.MOD(RefDet(z),2).eq.1) then
        iop=1
        Iopen(iop)=RefDet(z)
        ClosedShellCase=.false.
      endif
    enddo
    iMax=iop

    if(ClosedShellCase) then
      write(iout,*) 'Sorting coefficients in CLOSED-SHELL SYSTEM'
    else
      write(iout,*) 'Sorting coefficients in OPEN-SHELL SYSTEM'
      if (MOD(Iopen(iMax),2).eq.1) spo=1
      if (MOD(Iopen(iMax),2).eq.0) spo=-1
    endif


   ! open-shell systems require a specific reading in order to reorganize
   ! the alpha and beta spin orbitals to make it readable in Molpro as
   ! efficiently as possible. Here 2 different ways are implemented:
    if(sorting_way.eq.1) then

      write(iout,*) 'Coefficients listed in way 1:'
      write(iout,*) '  OCC(alpha),OCC(beta),VIR(alpha),VIR(beta)'

      ial=0
      ! loop to find all the alpha electrons
      do z=1,nel
        if (MOD(RefDet(z),2).eq.1) then
          ial=ial+1
          Itot(ial)=RefDet(z)
        endif
      enddo
      ialMax=ial
      ibe=0
      ! loop to find all the beta electrons
      do z=1,nel
        if (MOD(RefDet(z),2).eq.0) then
          ibe=ibe+1
          Itot(ialMax+ibe)=RefDet(z)
        endif
      enddo
      ibeMax=ibe
      if(ialMax+ibeMax.ne.nel) write(iout,*) 'WARNING: not matching number of electrons!'

      ialVir=0
      ! loop to find all the non-occupied alpha spin orbitals
      do i=1,nbasis
        j=0
        do z=1,nel
          if(i.eq.RefDet(z)) j=2
        enddo
        if(j.eq.2) cycle
        if (MOD(i,2).eq.1) then
          ialVir=ialVir+1
          Itot(nel+ialVir)=i
        endif
      enddo
      ialVirMax=ialVir
      ibeVir=0
      ! loop to find all the non-occupied beta spin orbitals
      do i=1,nbasis
        j=0
        do z=1,nel
          if(i.eq.RefDet(z)) j=2
        enddo
        if(j.eq.2) cycle
        if (MOD(i,2).eq.0) then
          ibeVir=ibeVir+1
          Itot(nel+ialVirMax+ibeVir)=i
        endif
      enddo
      ibeVirMax=ibeVir
      if(nel+ialVirMax+ibeVirMax.ne.nbasis) write(iout,*) 'WARNING: not matching number of orbitals!'

    else if(sorting_way.eq.2) then

      write(iout,*) 'Coefficients listed in way 2:'
      if(ClosedShellCase) then
        write(iout,*) '  CLOSED(alpha,beta),VIRTUALS(alpha,beta)'
      else
        write(iout,*) '  CLOSED(alpha,beta),OPEN(alpha),OPEN(beta),VIRTUALS(alpha,beta)'
      endif
    endif


      ! LOOP TO READ ALL THE CI COEFFICIENTS FROM ASCII FILES
      do Nind=1,n_store_ci_level
        write(fileCICoeffAv, '( "CI_COEFF_",I1,"_AV" )') Nind
        open (unit=30+Nind,file=fileCICoeffAv,status='old', action='read')
        write(fileCIcoeffSort, '("CI_COEFF_",I1)') Nind
        open (unit=100+Nind,file=fileCIcoeffSort,status='replace')
        do
          read(30+Nind,*,IOSTAT=z) x,(indCoef(1,i),indCoef(2,i),i=1,Nind)
          if (z<0) then
             exit
          endif
          signCI=1

          if(sorting_way.eq.1) then

            do i=1,Nind
               do p=1,nel
                  if(indCoef(1,i).eq.Itot(p)) then
                     indCoef(1,i)=p
                     exit
                  endif
               enddo
               do p=nel+1,nbasis
                  if(indCoef(2,i).eq.Itot(p)) then
                     indCoef(2,i)=p
                     exit
                  endif
               enddo
            enddo
          ! INSERTION SORT
            do p=1,2
              do i=2,Nind
                do j=i,2,-1
                  if (indCoef(p,j).lt.indCoef(p,j-1)) then
                    call swap(indCoef(p,j),indCoef(p,j-1))
                    signCI=-signCI
                  else
                    exit
                  endif
                enddo
              enddo
            enddo

          else if(sorting_way.eq.2) then
    ! TODO: in this way 2. of listing is not included yet the sorting of the open-shell
    !       RefDet with one or more virtual spatial orbs in between two occupied
    !       spin orbitals, e.g., p=4,6,8,...

            if(.not.ClosedShellCase) then

              do i=1,Nind
                do iop=1,iMax
                  if(indCoef(1,i).eq.Iopen(iop)) then
                    indCoef(1,i)=nel-iMax+iop
                    exit
                  endif
                enddo
                do iop=1,iMax
                  if(indCoef(2,i).eq.Iopen(iop)+spo) then
                    indCoef(2,i)=nel+iop
                    exit
                  endif
                enddo
              enddo
            endif
          endif

          if(Nind.eq.1) then
            S(indCoef(1,1),indCoef(2,1)) = x
          else if(Nind.eq.2) then
            D(indCoef(1,1),indCoef(2,1),indCoef(1,2),indCoef(2,2)) = signCI*x
          else if(Nind.eq.3) then
!            (T(indCoef(i,1),indCoef(i,2)),i=1,Nind) = signCI*x
            T(indCoef(1,1),indCoef(2,1),indCoef(1,2),indCoef(2,2),indCoef(1,3),indCoef(2,3)) = signCI*x
          endif
        enddo
        close (30+Nind)
      enddo


   ! Writing the CI coefficients in the ASCII files
    do i = 1,nel
      do a = nel+1,nbasis
        if(.not. near_zero(S(i,a))) then
          write(101,'(G20.12,2I5)') S(i,a),i,a
!          write(101,'(G20.12,2I5)') S(i,a),i,a-nel
        endif
      enddo
    enddo

    do j = 2,nel
      do i = 1, nel-1
        do a = nel+1,nbasis
          do b = a+1, nbasis
            if(.not. near_zero(D(i,a,j,b)).and. i.lt.j) then
              write(102,'(G20.12,4I5)') D(i,a,j,b),i,a,j,b
!              write(102,'(G20.12,4I5)') D(i,a,j,b),i,a-nel,j,b-nel
            endif
          enddo
        enddo
      enddo
    enddo

    do k = 3,nel
      do j = 2, nel-1
        do i = 1, nel-2
          do a = nel+1,nbasis
            do b = a+1, nbasis
              do c = b+1, nbasis
                if(.not. near_zero(T(i,a,j,b,k,c)).and. i.lt.j .and. j.lt.k) then
                  write(103,'(G20.12,6I5)') T(i,a,j,b,k,c),i,a,j,b,k,c
!                  write(103,'(G20.12,6I5)') T(i,a,j,b,k,c),i,a-nel,j,b-nel,k,c-nel
                end if
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo

    close(101)
    close(102)
    close(103)

  end subroutine sorting

  subroutine swap(a,b)
    integer,intent(inout) :: a,b
    integer                :: c
    c=a
    a=b
    b=c
  end subroutine swap

end module sdt_amplitudes
