! Module to collect and average the CI coefficients
module sdt_amplitudes

  use bit_reps, only: extract_sign, decode_bit_det, encode_sign, niftot, nifd
  use constants, only: dp, lenof_sign, n_int, int64, stdout
  use DetBitOps, only: get_bit_excitmat, EncodeBitDet, GetBitExcitation!, FindBitExcitLevel
  use util_mod, only: near_zero
  use FciMCData, only: TotWalkers, iLutRef, CurrentDets, AllNoatHf, projedet, &
                       ll_node
  use hash, only: hash_table_lookup, add_hash_table_entry, init_hash_table, &
                  clear_hash_table
  use LoggingData, only: n_store_ci_level, n_iter_ci_coeff
  use SystemData, only: nel, nbasis, symmax
  use Parallel_neci, only: iProcIndex, MPIcollection
  use MPI_wrapper, only: root

  implicit none

  integer(n_int), allocatable :: ciCoeff_storage(:,:), root_ciCoeff_storage(:,:)
  integer :: hash_table_ciCoeff_size, first_free_entry, nCyc, root_first_free_entry
  type(ll_node), pointer :: hash_table_ciCoeff(:)
  character (len=90) :: fileCICoeffAv, fileCIcoeffSort
  character (len=90) :: fileCICoeffTest, filecoefTestSrt0, filecoefTestSrt1
  character (len=90) :: filecoefTestSrt2, filecoefTestSrt3, filecoefTestSrt4
  character (len=90) :: filecoefTestSrt5, filecoefTestSrt6
  integer(n_int), allocatable  :: totex_coeff(:,:)

contains


  subroutine init_ciCoeff
    nCyc = 0
    first_free_entry = 0
    hash_table_ciCoeff_size = 500000
    allocate(hash_table_ciCoeff(hash_table_ciCoeff_size))
    call init_hash_table(hash_table_ciCoeff)
    allocate(ciCoeff_storage(0:NIfTot,hash_table_ciCoeff_size))
    allocate(totex_coeff(n_store_ci_level,n_store_ci_level))
  end subroutine init_ciCoeff


  subroutine storeCiCoeffs
    integer :: ic, ex(2,4),nIEx(nel)
    integer(int64) :: i
    real(dp) :: sign_tmp(lenof_sign)
    nCyc = nCyc + 1

    ! loop through all occupied determinants
    do i = 1, TotWalkers
       ! definition of the max ic as input for get_bit_excitmat
       ic = n_store_ci_level+1
!       ic = FindBitExcitLevel(ilutRef(:,1),CurrentDets(:,i))
       ! extraction of the excitation level from every determinant
       call get_bit_excitmat(iLutRef(:,1),CurrentDets(:,i),ex,ic)
       if(ic.le.n_store_ci_level) then
          call extract_sign(CurrentDets(:,i), sign_tmp)
          call decode_bit_det(nIEx, CurrentDets(:,i))
          call cache_sign(sign_tmp, nIEx)
       end if
    enddo

  end subroutine storeCiCoeffs


  ! it updates the CI coeffs storage list
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
!           write(stdout,*) 'sign_tmp = sign_tmp + sgn = ', sign_tmp, sgn
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


  ! it prints averaged CI coeffs collected during the calcualtion
  subroutine print_averaged_ci_coeff

    use SymExcitDataMod, only:  OrbClassCount
    use GenRandSymExcitNUMod , only : ClassCountInd
    use util_mod, only: intswap
    use sort_mod, only: sort

    integer :: i, ic, ex(2,4),icI, Nind,h1,h2,h3
    real(dp) :: sign_tmp(lenof_sign)
    logical  :: tPar
    integer  :: h,j,k,z,p,ial(symmax),ibe(symmax),Itot(nbasis),iSym,totEl,totOrb
    integer  :: signCI,indCoef(2,n_store_ci_level),Norb(symmax),NorbTot(symmax)

    if(iProcIndex.eq.root) then
      write(stdout,*) ''
      write(stdout,*) '-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --'
      write(stdout,*) ''
      write(stdout,*) '*** CI COEFFICIENTS COLLECTION ***'
      write(stdout,*) ''
      write(stdout,"(A44,I10)") 'Maximum excitation level of the CI coeffs = ', n_store_ci_level
      write(stdout,"(A44,I10)") 'Number of iterations set for average      = ', n_iter_ci_coeff
      if(nCyc.ne.n_iter_ci_coeff) then
         write(stdout,"(A44,I10)") ' -> actual iterations used for average    = ', nCyc
      endif
    endif

    call MPIcollection(NIfTot,first_free_entry,ciCoeff_storage,root_first_free_entry,root_ciCoeff_storage)

    if(iProcIndex.eq.root) then
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! loop to find all the alpha/beta occ/unocc orbs !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      NorbTot(:)=0
      do iSym = 1, symmax
        Norb(iSym) =  OrbClassCount(ClassCountInd(1/2,iSym,0)) * 2
        NorbTot(iSym) = Norb(iSym) + NorbTot(iSym-1)
      enddo
     ! loop to find all the occupied alpha spin orbitals
      totEl = 0
      do iSym = 1, symmax
        ial(iSym)=0
        do z=1,nel
          if(projEDet(z,1).gt.NorbTot(iSym-1) .and. projEDet(z,1).le.NorbTot(iSym)) then
            if (MOD(projEDet(z,1),2).eq.1) then
              ial(iSym)=ial(iSym)+1
              Itot(totEl+ial(iSym))=projEDet(z,1)
            endif
          endif
        enddo
        totEl = totEl + ial(iSym)
        ibe(iSym)=0
     ! loop to find all the occupied beta spin orbitals
        do z=1,nel
          if(projEDet(z,1).gt.NorbTot(iSym-1) .and. projEDet(z,1).le.NorbTot(iSym)) then
            if (MOD(projEDet(z,1),2).eq.0) then
              ibe(iSym)=ibe(iSym)+1
              Itot(totEl+ibe(iSym))=projEDet(z,1)
            endif
          endif
        enddo
        totEl = totEl + ibe(iSym)
      enddo
      if(totEl.ne.nel) write(stdout,*) 'WARNING: not matching number of electrons!'
     ! loop to find all the non-occupied alpha spin orbitals
      totOrb=0
      do iSym = 1, symmax
        ial(iSym)=0
        do k=NorbTot(iSym-1)+1,NorbTot(iSym)
          j=0
          do z=1,nel
            if(k.eq.projEDet(z,1)) j=2
          enddo
          if(j.eq.2) cycle
          if (MOD(k,2).eq.1) then
            ial(iSym)=ial(iSym)+1
            Itot(nel+totOrb+ial(iSym))=k
          endif
        enddo
        totOrb = totOrb + ial(iSym)
        ibe(iSym)=0
     ! loop to find all the non-occupied beta spin orbitals
        do k=NorbTot(iSym-1)+1,NorbTot(iSym)
          j=0
          do z=1,nel
            if(k.eq.projEDet(z,1)) j=2
          enddo
          if(j.eq.2) cycle
          if (MOD(k,2).eq.0) then
            ibe(iSym)=ibe(iSym)+1
            Itot(nel+totOrb+ibe(iSym))=k
          endif
        enddo
        totOrb = totOrb + ibe(iSym)
      enddo
      if(nel+totOrb.ne.nbasis) write(stdout,*) 'WARNING: not matching number of orbitals!'
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      totex_coeff(:,:) = 0
      do icI = 0, n_store_ci_level
         if(icI.ne.0) then
           write(fileCICoeffAv, '( "ci_coeff_",I1,"_av" )') icI
           open (unit=30+icI,file=fileCICoeffAv,status='replace')
           write(fileCICoeffTest, '( "ci_coeff_",I1,"_test" )') icI
           open (unit=40+icI,file=fileCICoeffTest,status='replace')
         endif
         ! loop over the total entries of CI coefficients
         do i = 1, root_first_free_entry
            signCI=1
            ic = n_store_ci_level+1
            ! call to get the excitation level of the CI coefficient
            call get_bit_excitmat(iLutRef(:,1),root_ciCoeff_storage(:,i),ex,ic)
!            ic = FindBitExcitLevel(ilutRef(:,1),root_ciCoeff_storage(:,i))
            if(icI.eq.ic) then
              ! call to get the numerical value of the CI coefficient
              ! (i.e. the number of walkers not yet averaged)
              call extract_sign(root_ciCoeff_storage(:,i),sign_tmp)
              ex(1,1) = ic
              ! call to get the sign of the CI coefficient
              ! if tPar is true, there is an odd number of permutations
              call GetBitExcitation(ilutRef(:,1),root_ciCoeff_storage(:,i),ex,tPar)
              if(tPar) sign_tmp = -sign_tmp

              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              ! indices conversion
               do k=1,icI
                 do p=1,nel
                     if(ex(1,k).eq.Itot(p)) then
                       indCoef(1,k)=p
                       exit
                    endif
                 enddo
                 do p=nel+1,nbasis
                    if(ex(2,k).eq.Itot(p)) then
                       indCoef(2,k)=p-nel
                       exit
                     endif
                 enddo
               enddo
               ! insertion sort
               do p=1,2
                 do k=2,icI
                   do j=k,2,-1
                     if (indCoef(p,j).lt.indCoef(p,j-1)) then
                       call intswap(indCoef(p,j),indCoef(p,j-1))
                       signCI=-signCI
                     else
                       exit
                     endif
                   enddo
                 enddo
               enddo
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

               select case(icI)
               case(0)
                  write(stdout,"(A44,F14.3)") 'Instantaneous number of walkers on HF     = ', AllNoatHf
                  AllNoatHf = -sign_tmp/nCyc
                  write(stdout,"(A44,F14.3)") 'Averaged number of walkers on HF          = ', AllNoatHf
                  write(stdout,"(A44,I10)") 'Total entries of CI coefficients          = ', root_first_free_entry
!                  write(stdout,*) 'sign_tmp/(AllNoatHf*nCyc) =', sign_tmp/(AllNoatHf*nCyc)
               case(1)
                 totex_coeff(icI,1) = totex_coeff(icI,1) + 1
                 write(31,'(G20.12,2I5)') sign_tmp/(AllNoatHf(1)*nCyc), &
                                           ex(1,1),ex(2,1)
                 if(.not. near_zero(sign_tmp(1))) then
                   totex_coeff(icI,2) = totex_coeff(icI,2) + 1
                   write(41,'(G20.12,2I5)') signCI*(sign_tmp/(AllNoatHf(1)*nCyc)),&
                                            indCoef(1,1),indCoef(2,1)
                 endif
               case(2)
                 totex_coeff(icI,1) = totex_coeff(icI,1) + 1
                 write(32,'(G20.12,4I5)') sign_tmp/(AllNoatHf(1)*nCyc),ex(1,1),&
                                           ex(2,1),ex(1,2),ex(2,2)
                 if(.not. near_zero(sign_tmp(1))) then
                   totex_coeff(icI,2) = totex_coeff(icI,2) + 1
                   write(42,'(G20.12,4I5)') signCI*(sign_tmp(1)/(AllNoatHf(1)*nCyc)),&
                               indCoef(1,1),indCoef(2,1),indCoef(1,2),indCoef(2,2)
                 endif
               case(3)
                 totex_coeff(icI,1) = totex_coeff(icI,1) + 1
                 write(33,'(G20.12,6I5)') sign_tmp/(AllNoatHf(1)*nCyc),ex(1,1),&
                                           ex(2,1),ex(1,2),ex(2,2),ex(1,3),ex(2,3)
                 if(.not. near_zero(sign_tmp(1))) then
                   totex_coeff(icI,2) = totex_coeff(icI,2) + 1
                   write(43,'(G20.12,6I5)') signCI*(sign_tmp/(AllNoatHf(1)*nCyc)),&
            indCoef(1,1),indCoef(2,1),indCoef(1,2),indCoef(2,2),indCoef(1,3),indCoef(2,3)
                 endif
               end select
            end if
         enddo
         if(iCI.ne.0) close(30+icI)
         if(iCI.ne.0) close(40+icI)
      enddo
      write(stdout,"(A44,I10)") ' total entries for singles                = ', totex_coeff(1,1)
      if(totex_coeff(1,1).ne.totex_coeff(1,2)) then
        write(stdout,"(A44,I10)") ' -total entries for singles without zeros = ',totex_coeff(1,2)
      endif
      write(stdout,"(A44,I10)") ' total entries for doubles                = ', totex_coeff(2,1)
      if(totex_coeff(2,1).ne.totex_coeff(2,2)) then
        write(stdout,"(A44,I10)") ' -total entries for doubles without zeros = ',totex_coeff(2,2)
      endif
      if(n_store_ci_level.eq.3) then
        write(stdout,"(A44,I10)") ' total entries for triples                = ', totex_coeff(3,1)
        if(totex_coeff(3,1).ne.totex_coeff(3,2)) then
          write(stdout,"(A44,I10)") ' -total entries for triples without zeros = ',totex_coeff(3,2)
        endif
      endif

      write(stdout,*) ' sorting CI coefficients...'
      call sorting()

      write(stdout,*) '-> CI coefficients written in ASCII files ci_coeff_*'
      write(stdout,*) ''
      write(stdout,*) '-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --'
    endif

    call fin_ciCoeff()
  end subroutine print_averaged_ci_coeff


  subroutine fin_ciCoeff
    call clear_hash_table(hash_table_ciCoeff)
    deallocate(hash_table_ciCoeff)
    deallocate(ciCoeff_storage)
  end subroutine fin_ciCoeff


  ! it lists the averaged CI coeffs sorting the indices in
  ! this way: OCC(alpha), OCC(beta), VIR(alpha), VIR(beta)
  subroutine sorting

  type :: t_singles
    integer :: a,i
    double precision :: x
  end type
  type :: t_doubles
    integer :: a,b,i,j
    double precision :: x
  end type
  type :: t_triples
    integer :: a,b,c,i,j,k
    double precision :: x
  end type

  integer :: hI,Nind,z
  integer :: a,b,c,j,i,k
  type(t_singles),allocatable :: singles(:), singles_tmp(:)
  type(t_doubles),allocatable :: doubles(:), doubles_tmp(:)
  type(t_triples),allocatable :: triples(:), triples_tmp(:)

  do Nind=1,n_store_ci_level
    write(fileCICoeffTest, '( "ci_coeff_",I1,"_test" )') Nind
    open (unit=40+Nind,file=fileCICoeffTest,status='old', action='read')
    write(filecoefTestSrt0, '("ci_coeff_",I1 )') Nind
    open (unit=140+Nind,file=filecoefTestSrt0,status='replace')

    if(Nind.eq.1) allocate(singles(totex_coeff(1,2)))
    if(Nind.eq.2) allocate(doubles(totex_coeff(2,2)))
    if(Nind.eq.3) allocate(triples(totex_coeff(3,2)))
    if(Nind.eq.1) allocate(singles_tmp(totex_coeff(1,2)))
    if(Nind.eq.2) allocate(doubles_tmp(totex_coeff(2,2)))
    if(Nind.eq.3) allocate(triples_tmp(totex_coeff(3,2)))

    hI=0
    do
     hI = hI + 1
     if(Nind.eq.1) then
       read(40+Nind,*,IOSTAT=z) singles(hI)%x,singles(hI)%i,&
                                singles(hI)%a
     else if(Nind.eq.2) then
       read(40+Nind,*,IOSTAT=z) doubles(hI)%x,doubles(hI)%i,&
                  doubles(hI)%a,doubles(hI)%j,doubles(hI)%b
     else if(Nind.eq.3) then
       read(40+Nind,*,IOSTAT=z) triples(hI)%x,triples(hI)%i,&
                  triples(hI)%a,triples(hI)%j,triples(hI)%b,&
                                triples(hI)%k,triples(hI)%c
     endif
     if (z<0) exit
    enddo
    close (40+Nind, status='delete')


   ! sorting singles
   if(Nind.eq.1) then
    z=0
    do a = 1,nbasis-nel
      do hI = 1,totex_coeff(Nind,2)
       if(singles(hI)%a.eq.a) then
        z = z + 1
        singles_tmp(z) = singles(hI)
       endif
      enddo
    enddo
    do i = 1,nel
      do hI = 1,totex_coeff(Nind,2)
       if(singles_tmp(hI)%i.eq.i) then
        write(140+Nind,'(G20.12,2I5)') singles_tmp(hI)%x,singles_tmp(hI)%i,&
                                       singles_tmp(hI)%a
       endif
      enddo
    enddo
   close (140+Nind)


   ! sorting doubles
   else if(Nind.eq.2) then
    z=0
    do a = 1,nbasis-nel-1
     do hI = 1,totex_coeff(Nind,2)
      if(doubles(hI)%a.eq.a) then
        z = z + 1
        doubles_tmp(z) = doubles(hI)
       endif
     enddo
    enddo
    z=0
    do b = 2,nbasis-nel
      do hI = 1,totex_coeff(Nind,2)
      if(doubles_tmp(hI)%b.eq.b) then
        z = z + 1
        doubles(z) = doubles_tmp(hI)
        endif
      enddo
    enddo
    z=0
    do i = 1,nel-1
      do hI = 1,totex_coeff(Nind,2)
      if(doubles(hI)%i.eq.i) then
        z = z + 1
        doubles_tmp(z) = doubles(hI)
        endif
      enddo
    enddo
    do j = 2,nel
      do hI = 1,totex_coeff(Nind,2)
      if(doubles_tmp(hI)%j.eq.j) then
       write(140+Nind,'(G20.12,4I5)') doubles_tmp(hI)%x,doubles_tmp(hI)%i,&
                    doubles_tmp(hI)%a,doubles_tmp(hI)%j,doubles_tmp(hI)%b
        endif
      enddo
    enddo
   close (140+Nind)


   ! sorting triples
   else if(Nind.eq.3) then
    z=0
    do c = 3,nbasis-nel
      do hI = 1,totex_coeff(Nind,2)
      if(triples(hI)%c.eq.c) then
       z = z + 1
       triples_tmp(z) = triples(hI)
      endif
     enddo
    enddo
    z=0
    do b = 2,nbasis-nel-1
      do hI = 1,totex_coeff(Nind,2)
      if(triples_tmp(hI)%b.eq.b) then
       z=z+1
       triples(z) = triples_tmp(hI)
      endif
     enddo
    enddo
    z=0
    do a = 1,nbasis-nel-2
      do hI = 1,totex_coeff(Nind,2)
      if(triples(hI)%a.eq.a) then
       z=z+1
       triples_tmp(z) = triples(hI)
      endif
     enddo
    enddo
    z=0
    do i = 1,nel-2
      do hI = 1,totex_coeff(Nind,2)
      if(triples_tmp(hI)%i.eq.i) then
       z=z+1
       triples(z) = triples_tmp(hI)
       endif
      enddo
    enddo
    z=0
    do j = 2,nel-1
      do hI = 1,totex_coeff(Nind,2)
      if(triples(hI)%j.eq.j) then
       z=z+1
       triples_tmp(z) = triples(hI)
      endif
     enddo
    enddo
    do k = 3,nel
      do hI = 1,totex_coeff(Nind,2)
      if(triples_tmp(hI)%k.eq.k) then
       write(140+Nind,'(G20.12,6I5)') triples_tmp(hI)%x,triples_tmp(hI)%i,&
                    triples_tmp(hI)%a,triples_tmp(hI)%j,triples_tmp(hI)%b,&
                                      triples_tmp(hI)%k,triples_tmp(hI)%c
      endif
     enddo
    enddo
   close (140+Nind)
   endif

  if(Nind.eq.1) deallocate(singles)
  if(Nind.eq.2) deallocate(doubles)
  if(Nind.eq.3) deallocate(triples)
  if(Nind.eq.1) deallocate(singles_tmp)
  if(Nind.eq.2) deallocate(doubles_tmp)
  if(Nind.eq.3) deallocate(triples_tmp)

  enddo ! do Nind=1,n_store_ci_level

  end subroutine sorting

end module sdt_amplitudes
