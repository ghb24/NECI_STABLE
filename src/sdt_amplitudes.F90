module sdt_amplitudes

  use bit_reps, only: extract_sign, decode_bit_det, encode_sign, niftot, nifd
  use constants, only: dp, lenof_sign, EPS, n_int, bits_n_int
  use DetBitOps, only: get_bit_excitmat, FindBitExcitLevel, EncodeBitDet
  use FciMCData, only: TotWalkers, iLutRef, CurrentDets, AllNoatHf, projedet, &
                       HashIndex, CurrentDets, ll_node, norm_psi, ilutHF, &
                       VaryShiftIter, iter
  use hash, only: hash_table_lookup,add_hash_table_entry,init_hash_table, &
                  clear_hash_table
  use LoggingData, only: n_store_ci_level
  use SystemData, only: nel,nbasis
  
  implicit none

  integer(n_int), allocatable :: ciCoeff_storage(:,:)
  integer :: hash_table_ciCoeff_size, first_free_entry = 0
  type(ll_node), pointer :: hash_table_ciCoeff(:)
  integer :: iout = 6
  integer :: nCyc,nREF,nS,nD,nT, storSize
  real(dp) :: ciCoeff_REF(lenof_sign)
  real(dp), allocatable :: ciCoeff_storage_S(:,:)
  real(dp), allocatable :: ciCoeff_storage_D(:,:,:,:)
  real(dp), allocatable :: ciCoeff_storage_T(:,:,:,:,:,:)

contains


  subroutine init_ciCoeff
    nCyc = 0
    first_free_entry = 0
    hash_table_ciCoeff_size = 50000
    storSize = 20
    call init_hash_table(hash_table_ciCoeff)
    allocate(hash_table_ciCoeff(hash_table_ciCoeff_size))
    allocate(ciCoeff_storage(0:NIfTot,hash_table_ciCoeff_size))
  end subroutine init_ciCoeff


  subroutine init_storeCiCoeff
    allocate(ciCoeff_storage_S(storSize,storSize))
    allocate(ciCoeff_storage_D(storSize,storSize,storSize,storSize))
    allocate(ciCoeff_storage_T(storSize,storSize,storSize,storSize,storSize,storSize))
    ciCoeff_REF = 0.d0
    ciCoeff_storage_S = 0.d0
    ciCoeff_storage_D = 0.d0
    ciCoeff_storage_T = 0.d0
  end subroutine init_storeCiCoeff


  subroutine fin_ciCoeff
    call clear_hash_table(hash_table_ciCoeff)
    deallocate(hash_table_ciCoeff)
    deallocate(ciCoeff_storage)
  end subroutine fin_ciCoeff


  subroutine print_snapshot_ci_coeff
    integer :: i, ic, ex(2,3),icI
    real(dp) :: sign_tmp(lenof_sign)

    open (unit=21,file='SINGLES',status='replace')
    open (unit=22,file='DOUBLES',status='replace')
    open (unit=23,file='TRIPLES',status='replace')
    write(21,*) '&FCI NBASIS=', nbasis,',', 'NELEC=', nel
    write(21,*) '&END'
    write(22,*) '&FCI NBASIS=', nbasis,',', 'NELEC=', nel
    write(22,*) '&END'
    write(23,*) '&FCI NBASIS=', nbasis,',', 'NELEC=', nel
    write(23,*) '&END'

    write(iout,*) 'CI level =', n_store_ci_level

    do icI = 0, n_store_ci_level
      do i = 1, TotWalkers
        ic = 4
        call get_bit_excitmat(iLutRef(:,1), CurrentDets(:,i),ex,ic)
        call extract_sign(CurrentDets(:,i), sign_tmp)
        if(icI.eq.ic) then
          select case(icI)
          case(0)
            AllNoatHf = sign_tmp
            write(iout,*) 'REFERENCE:      - N el.=',nel,'N basis=',nbasis
            write(iout,*) 'Total Walkers = ', TotWalkers
            write(iout,*) sign_tmp,AllNoatHf,sign_tmp(1),AllNoatHf(1), &
                          sign_tmp(1)/AllNoatHf(1)
          case(1)
!!            write(21,'(A10,G20.12,2I5)') 'SINGLES: ', sign_tmp/AllNoatHf(1),ex(1,1),ex(2,1)
            write(21,'(G20.12,2I5)') sign_tmp/AllNoatHf(1),ex(1,1),ex(2,1)
          case(2)
            write(22,'(G20.12,4I5)') sign_tmp/AllNoatHf(1),ex(1,1),ex(2,1),&
                                     ex(1,2),ex(2,2)
          case(3)
            write(23,'(G20.12,6I5)') sign_tmp/AllNoatHf(1),ex(1,1),ex(2,1),&
                                     ex(1,2),ex(2,2),ex(1,3),ex(2,3)
          end select
        end if
      end do
    end do

    close(21)
    close(22)
    close(23)
  end subroutine print_snapshot_ci_coeff


  subroutine print_averaged_ci_coeff
    integer :: i, ic, ex(2,3), icI
    real(dp) :: sign_tmp(lenof_sign)

!    to sort the excitations i need to
!    call sort(ciCoeff_storage)
!    present in quickLib

    open (unit=31,file='SINGLES-AV',status='replace')
    open (unit=32,file='DOUBLES-AV',status='replace')
    open (unit=33,file='TRIPLES-AV',status='replace')
!    write(31,*) '&FCI NBASIS=', nbasis,',', 'NELEC=', nel
!    write(31,*) '&END'
!    write(32,*) '&FCI NBASIS=', nbasis,',', 'NELEC=', nel
!    write(32,*) '&END'
!    write(33,*) '&FCI NBASIS=', nbasis,',', 'NELEC=', nel
!    write(33,*) '&END'

    do icI = 0, n_store_ci_level
       do i = 1, first_free_entry
          ic = 4
          call get_bit_excitmat(iLutRef(:,1),ciCoeff_storage(:,i),ex,ic)
          if(icI.eq.ic) then
             select case(icI)
             case(0)
                call extract_sign(ciCoeff_storage(:,i),sign_tmp)
                AllNoatHf = sign_tmp/nCyc
                write(iout,*) 'REFERENCE-AV:    - N el.=',nel,'N basis=',nbasis
                write(iout,*) 'Total S+D+T = ', first_free_entry, 'nCyc=', nCyc
                write(iout,*) sign_tmp/nCyc,AllNoatHf,sign_tmp(1),AllNoatHf(1),&
                     (sign_tmp(1)/nCyc)/AllNoatHf(1)
             case(1)
                call extract_sign(ciCoeff_storage(:,i),sign_tmp)
                write(31,'(G20.12,2I5)') sign_tmp/(AllNoatHf(1)*nCyc), & 
                                         ex(1,1),ex(2,1)
             case(2) 
                call extract_sign(ciCoeff_storage(:,i),sign_tmp)
!                write(iout,'(A)',advance='no') 'DOUBLES-AV:'
                write(32,'(G20.12,4I5)') sign_tmp/(AllNoatHf(1)*nCyc),ex(1,1),&
                                         ex(2,1),ex(1,2),ex(2,2)
             case(3)
                call extract_sign(ciCoeff_storage(:,i),sign_tmp)
!                write(iout,'(A)',advance='no') 'TRIPLES-AV:'
                write(33,'(G20.12,6I5)') sign_tmp/(AllNoatHf(1)*nCyc),ex(1,1),&
                                         ex(2,1),ex(1,2),ex(2,2),ex(1,3),ex(2,3)
             end select
          end if
       end do
    end do

    close(31)
    close(32)
    close(33)
    call sorting()
    call fin_ciCoeff()
  end subroutine print_averaged_ci_coeff


  subroutine print_storeCiCoefficients
    integer :: i,j,k,l,a,b,c, ic, ex(2,3),REF
    real(dp) :: sign_tmp(lenof_sign)

    open (unit=31,file='SINGLES-AV',status='replace')
    open (unit=32,file='DOUBLES-AV',status='replace')
    open (unit=33,file='TRIPLES-AV',status='replace')
    write(31,*) '&FCI NBASIS=', nbasis,',', 'NELEC=', nel
    write(31,*) '&END'
    write(32,*) '&FCI NBASIS=', nbasis,',', 'NELEC=', nel
    write(32,*) '&END'
    write(33,*) '&FCI NBASIS=', nbasis,',', 'NELEC=', nel
    write(33,*) '&END'

    write(iout,*) 'N el.=',nel,'N basis=',nbasis
    write(iout,*) 'Total Walkers = ', TotWalkers, 'nREF=', nREF
    write(iout,*) 'REFERENCE-AV: ',  ciCoeff_REF(1)/nREF
    write(iout,*) 'nSingles = ',  nS, 'nDoubles = ', nD, 'nTriples = ', nT
!    REF = ciCoeff_REF(1)/nREF
    REF = AllNoatHf(1)
    do i = 1,nel
      do a = nel+1,nbasis
       if(abs(ciCoeff_storage_S(i,a)).lt.0.00000001) then 
         ciCoeff_storage_S(i,a)=0.d0
       end if
!!       write(31,'(A10,G20.12,2I5)') 'SING-AV: ', (ciCoeff_storage_S(i,a)/nREF)/REF,i,a
       write(31,'(G20.12,2I5)') (ciCoeff_storage_S(i,a)/nREF)/REF,i,a
        do j = 1, nel
          do b = nel+1, nbasis
            if(abs(ciCoeff_storage_D(i,a,j,b)).lt.0.00000001) then 
              ciCoeff_storage_D(i,a,j,b)=0.d0
            end if
            write(32,'(G20.12,4I5)') (ciCoeff_storage_D(i,a,j,b)/nREF)/REF, & 
                                                      i,a,j,b
            do k = 1, nel
              do c = nel+1, nbasis
                if(abs(ciCoeff_storage_T(i,a,j,b,k,c)).lt.0.00000001) then 
                  ciCoeff_storage_T(i,a,j,b,k,c)=0.d0
                end if
              write(33,'(G20.12,6I5)') (ciCoeff_storage_T(i,a,j,b,k,c)/nREF)/REF, & 
                                       i,a,j,b,k,c
              end do
            end do
          enddo
        enddo
      enddo
    enddo

    close(31)
    close(32)
    close(33)
    call fin_ciCoeff()
  end subroutine print_storeCiCoefficients


  subroutine storeCiCoefficients
    integer :: i,j,k,l,a,b,c, ic, ex(2,3),nIEx(nel)
    real(dp) :: sign_tmp(lenof_sign)

    ! loop through all occ. determinants 
    do l = 1, TotWalkers
       ic = 4
       call get_bit_excitmat(iLutRef(:,1),CurrentDets(:,l),ex,ic)
       if(ic.eq.0) then
          call extract_sign(CurrentDets(:,l), sign_tmp)
          if(ciCoeff_REF(1).eq.0) then
            nREF=1
            ciCoeff_REF = sign_tmp
          else if(ciCoeff_REF(1).ne.0) then
            nREF=nREF+1
            ciCoeff_REF = ciCoeff_REF + sign_tmp
          end if
       else if(ic.eq.1) then
          call extract_sign(CurrentDets(:,l), sign_tmp)
          i = ex(1,1)
          a = ex(2,1)
          if(ciCoeff_storage_S(i,a).eq.0) then
            nS=1
            ciCoeff_storage_S(i,a) = sign_tmp(1)
          else if(ciCoeff_storage_S(i,a).ne.0) then
            nS=nS+1
            ciCoeff_storage_S(i,a) = ciCoeff_storage_S(i,a) + sign_tmp(1)
          end if
       else if(ic.eq.2) then
          call extract_sign(CurrentDets(:,l), sign_tmp)
          i = ex(1,1)
          a = ex(2,1)
          j = ex(1,2)
          b = ex(2,2)
          if(ciCoeff_storage_D(i,a,j,b).eq.0) then
            nD=1
            ciCoeff_storage_D(i,a,j,b) = sign_tmp(1)
          else if(ciCoeff_storage_D(i,a,j,b).ne.0) then
            call extract_sign(CurrentDets(:,l), sign_tmp)
            nD=nD+1
            ciCoeff_storage_D(i,a,j,b) = ciCoeff_storage_D(i,a,j,b) + sign_tmp(1)
          end if
       else if(ic.eq.3) then
          call extract_sign(CurrentDets(:,l), sign_tmp)
          i = ex(1,1)
          a = ex(2,1)
          j = ex(1,2)
          b = ex(2,2)
          k = ex(1,3)
          c = ex(2,3)
          if(ciCoeff_storage_T(i,a,j,b,k,c).eq.0) then
            nT=1
            ciCoeff_storage_T(i,a,j,b,k,c) = sign_tmp(1)
          else if((ic.eq.3).and.(ciCoeff_storage_T(i,a,j,b,k,c).ne.0)) then
            nT=nT+1
            ciCoeff_storage_T(i,a,j,b,k,c) = ciCoeff_storage_T(i,a,j,b,k,c) + sign_tmp(1)
          end if
       end if
    end do

  end subroutine storeCiCoefficients

  
  subroutine storeCiCoeffs
    integer :: i, ic, ex(2,3),nIEx(nel)
    real(dp) :: sign_tmp(lenof_sign)
    nCyc = nCyc + 1

    ! loop through all occupied determinants 
    do i = 1, TotWalkers
       ic = 4
       call get_bit_excitmat(iLutRef(:,1),CurrentDets(:,i),ex,ic)
       if(ic < 4) then
          call extract_sign(CurrentDets(:,i), sign_tmp)
          call decode_bit_det(nIEx, CurrentDets(:,i))
          call cache_sign(sign_tmp,nIEx,ex)
       end if
    end do

  end subroutine storeCiCoeffs



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
    if(tSuccess) then
       call extract_sign(ciCoeff_storage(:,ind),sign_tmp)
       sign_tmp = sign_tmp + sgn
!        if(all(nIEx.eq.projEDet(:,1))) then  !it is true only with the reference determinant
!           write(iout,*) 'sign_tmp = sign_tmp + sgn = ', sign_tmp, sgn
!        endif
    call encode_sign(ciCoeff_storage(:,ind),sign_tmp)
    ! else: add a new entry to our CI coeff storage
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


  subroutine sorting

    ! CAMBIARE SORTING COL NUOVO in READER2222.F90

    Implicit none
 
    double precision :: x
    double precision,allocatable :: S(:,:),D(:,:,:,:),T(:,:,:,:,:,:)
    integer :: i,j,k,a,b,c,z
 
    open (unit=101,file='SINGLES-AV',status='old', action='read')
    open (unit=102,file='DOUBLES-AV',status='old', action='read')
    open (unit=103,file='TRIPLES-AV',status='old', action='read')
    open (unit=31,file='SINGLES-AV-OR',status='replace')
    open (unit=32,file='DOUBLES-AV-OR',status='replace')
    open (unit=33,file='TRIPLES-AV-OR',status='replace')
 
!    write(31,*) '&FCI NBASIS=', nbasis,',', 'NELEC=', nel 
!    write(31,*) '&END'
!    write(32,*) '&FCI NBASIS=', nbasis,',', 'NELEC=', nel 
!    write(32,*) '&END'
!    write(33,*) '&FCI NBASIS=', nbasis,',', 'NELEC=', nel 
!    write(33,*) '&END'
 
    allocate(S(nel,nbasis))
    allocate(D(nel,nbasis,nel,nbasis))
    allocate(T(nel,nbasis,nel,nbasis,nel,nbasis))
    S(:,:)=0.d0
    D(:,:,:,:)=0.d0
    T(:,:,:,:,:,:)=0.d0
 
    do
      read(101,*,IOSTAT=z) x,i,a
      S(i,a) = x
     if (z<0) then
       exit
     endif
    enddo
    close (101)
 
    do
      read(102,*,IOSTAT=z) x,i,a,j,b
      D(i,a,j,b) = x
     if (z<0) then
       exit
     endif
    enddo
    close (102)
 
    do
      read(103,*,IOSTAT=z) x,i,a,j,b,k,c
      T(i,a,j,b,k,c) = x
     if (z<0) then
       exit
     endif
    enddo
    close (103)

    do i = 1,nel
      do a = nel+1,nbasis 
        if(abs(S(i,a)).ne.0) then 
          write(31,'(G20.12,2I5)') S(i,a),i,a-nel
        end if
      end do
    end do

    do j = 2,nel
      do i = 1, nel-1
        do a = nel+1,nbasis 
          do b = a+1, nbasis
            if(abs(D(i,a,j,b)).ne.0 .and. i.lt.j) then 
              write(32,'(G20.12,4I5)') D(i,a,j,b),i,a-nel,j,b-nel
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
                  write(33,'(G20.12,6I5)') T(i,a,j,b,k,c),i,a-nel,j,b-nel,k,c-nel
                end if
              enddo
            enddo
          enddo
        enddo
      enddo
    end do
 
  end subroutine sorting

end module sdt_amplitudes
