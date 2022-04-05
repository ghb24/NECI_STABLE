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
!  use sdt_util, only: sorting_singles

  implicit none

  integer(n_int), allocatable :: ciCoeff_storage(:,:), root_ciCoeff_storage(:,:)
  integer :: hash_table_ciCoeff_size, first_free_entry, nCyc, root_first_free_entry
  type(ll_node), pointer :: hash_table_ciCoeff(:)
  character (len=90) :: fileCICoeffAv, fileCIcoeffSort
  character (len=90) :: fileCICoeffTest, filecoefTestSrt0, filecoefTestSrt1
  character (len=90) :: filecoefTestSrt2, filecoefTestSrt3, filecoefTestSrt4
  character (len=90) :: filecoefTestSrt5, filecoefTestSrt6
  integer(n_int), allocatable  :: totex_coeff(:)

contains


  subroutine init_ciCoeff
    nCyc = 0
    first_free_entry = 0
    hash_table_ciCoeff_size = 500000
    allocate(hash_table_ciCoeff(hash_table_ciCoeff_size))
    call init_hash_table(hash_table_ciCoeff)
    allocate(ciCoeff_storage(0:NIfTot,hash_table_ciCoeff_size))
    allocate(totex_coeff(n_store_ci_level))
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
!    integer  :: totex_coeff(n_store_ci_level)

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
              NorbTot(0)=0
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

      totex_coeff(:) = 0
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
                 totex_coeff(icI) = totex_coeff(icI) + 1
                 write(31,'(G20.12,2I5)') sign_tmp/(AllNoatHf(1)*nCyc), &
                                           ex(1,1),ex(2,1)
                 write(41,'(G20.12,2I5)') signCI*(sign_tmp/(AllNoatHf(1)*nCyc)),&
                                            indCoef(1,1),indCoef(2,1)
               case(2)
                 totex_coeff(icI) = totex_coeff(icI) + 1
                 write(32,'(G20.12,4I5)') sign_tmp/(AllNoatHf(1)*nCyc),ex(1,1),&
                                           ex(2,1),ex(1,2),ex(2,2)
                 write(42,'(G20.12,4I5)') signCI*(sign_tmp/(AllNoatHf(1)*nCyc)),&
                             indCoef(1,1),indCoef(2,1),indCoef(1,2),indCoef(2,2)
               case(3)
                 totex_coeff(icI) = totex_coeff(icI) + 1
                 write(33,'(G20.12,6I5)') sign_tmp/(AllNoatHf(1)*nCyc),ex(1,1),&
                                           ex(2,1),ex(1,2),ex(2,2),ex(1,3),ex(2,3)
                 write(43,'(G20.12,6I5)') signCI*(sign_tmp/(AllNoatHf(1)*nCyc)),&
            indCoef(1,1),indCoef(2,1),indCoef(1,2),indCoef(2,2),indCoef(1,3),indCoef(2,3)
               end select
            end if
         enddo
         if(iCI.ne.0) close(30+icI)
         if(iCI.ne.0) close(40+icI)
      enddo
      write(stdout,"(A44,I10)") ' total entries for singles                = ', totex_coeff(1)
      write(stdout,"(A44,I10)") ' total entries for doubles                = ', totex_coeff(2)
      if(n_store_ci_level.eq.3) then
         write(stdout,"(A44,I10)") ' total entries for triples                = ', totex_coeff(3)
      endif

      write(stdout,*) ' sorting CI coefficients...'
      call sorting()
!      call sorting_singles(coeffx1,totex_coeff)

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! sorting by reading/writing !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
  type(t_singles),allocatable :: singles(:)
  type(t_doubles),allocatable :: doubles(:)
  type(t_triples),allocatable :: triples(:)

  allocate(singles(totex_coeff(1)))
  allocate(doubles(totex_coeff(2)))
  if(n_store_ci_level.eq.3) allocate(triples(totex_coeff(3)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do Nind=1,n_store_ci_level
    write(fileCICoeffTest, '( "ci_coeff_",I1,"_test" )') Nind
    open (unit=40+Nind,file=fileCICoeffTest,status='old', action='read')
    write(filecoefTestSrt0, '("ci_coeff_",I1,"_sort0" )') Nind
    open (unit=140+Nind,file=filecoefTestSrt0,status='replace')
    write(filecoefTestSrt1, '("ci_coeff_",I1,"_sort1" )') Nind
    open (unit=150+Nind,file=filecoefTestSrt1,status='replace')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! sorting by reading/writing from text files !!!!!!!!!!!!!!
    write(filecoefTestSrt2, '("ci_coeff_",I1,"_sort2" )') Nind
    open (unit=160+Nind,file=filecoefTestSrt2,status='replace')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if(Nind.ge.2) then
    write(filecoefTestSrt3, '("ci_coeff_",I1,"_sort3" )') Nind
    open (unit=170+Nind,file=filecoefTestSrt3,status='replace')
    write(filecoefTestSrt4, '("ci_coeff_",I1,"_sort4" )') Nind
    open (unit=180+Nind,file=filecoefTestSrt4,status='replace')
   endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if(Nind.eq.3) then
    write(filecoefTestSrt5, '("ci_coeff_",I1,"_sort5" )') Nind
    open (unit=190+Nind,file=filecoefTestSrt5,status='replace')
    write(filecoefTestSrt6, '("ci_coeff_",I1,"_sort6" )') Nind
    open (unit=200+Nind,file=filecoefTestSrt6,status='replace')
   endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    close (40+Nind)


    ! writing CI coefficients
    do hI = 1,totex_coeff(Nind)
     if(Nind.eq.1) then
      write(140+Nind,'(G20.12,2I5)') singles(hI)%x,singles(hI)%i,&
                                     singles(hI)%a
     else if(Nind.eq.2) then
      write(140+Nind,'(G20.12,4I5)') doubles(hI)%x,doubles(hI)%i,&
                       doubles(hI)%a,doubles(hI)%j,doubles(hI)%b
     else if(Nind.eq.3) then
      write(140+Nind,'(G20.12,6I5)') triples(hI)%x,triples(hI)%i,&
                       triples(hI)%a,triples(hI)%j,triples(hI)%b,&
                                     triples(hI)%k,triples(hI)%c
     endif
    enddo
    close (140+Nind)


   ! sorting singles
   if(Nind.eq.1) then
    do a = 1,nbasis-nel
      do hI = 1,totex_coeff(1)
       if(singles(hI)%a.eq.a.and..not.near_zero(singles(hI)%x)) then
        write(150+Nind,'(G20.12,2I5)') singles(hI)%x,singles(hI)%i,&
                                       singles(hI)%a
        endif
      enddo
    enddo
    close(150+Nind)

    open (unit=150+Nind,file=filecoefTestSrt1,status='old', action='read')
    hI=0
    do
      hI = hI + 1
      read(150+Nind,*,IOSTAT=z) singles(hI)%x,singles(hI)%i,&
                                singles(hI)%a
      if (z<0) exit
    enddo
    close(150+Nind)
    do i = 1,nel
      do hI = 1,totex_coeff(1)
       if(singles(hI)%i.eq.i) then
        write(160+Nind,'(G20.12,2I5)') singles(hI)%x,singles(hI)%i,&
                                       singles(hI)%a
       endif
      enddo
    enddo


   ! sorting doubles
   else if(Nind.eq.2) then
    do a = 1,nbasis-nel-1
     do hI = 1,totex_coeff(2)
       if(doubles(hI)%a.eq.a.and..not.near_zero(doubles(hI)%x)) then
        write(150+Nind,'(G20.12,4I5)') doubles(hI)%x,doubles(hI)%i,&
                         doubles(hI)%a,doubles(hI)%j,doubles(hI)%b
        endif
      enddo
    enddo
    close(150+Nind)

    open (unit=150+Nind,file=filecoefTestSrt1,status='old', action='read')
    hI=0
    do
      hI = hI + 1
      read(150+Nind,*,IOSTAT=z) doubles(hI)%x,doubles(hI)%i,&
                  doubles(hI)%a,doubles(hI)%j,doubles(hI)%b
      if (z<0) exit
    enddo
    close(150+Nind)
    do b = 2,nbasis-nel
     do hI = 1,totex_coeff(2)
      if(doubles(hI)%b.eq.b) then
       write(160+Nind,'(G20.12,4I5)') doubles(hI)%x,doubles(hI)%i,&
                        doubles(hI)%a,doubles(hI)%j,doubles(hI)%b
        endif
      enddo
    enddo
    close(160+Nind)

    open (unit=160+Nind,file=filecoefTestSrt2,status='old', action='read')
    hI=0
    do
      hI = hI + 1
      read(160+Nind,*,IOSTAT=z) doubles(hI)%x,doubles(hI)%i,&
                  doubles(hI)%a,doubles(hI)%j,doubles(hI)%b
      if (z<0) exit
    enddo
    close(160+Nind)
    do i = 1,nel-1
     do hI = 1,totex_coeff(2)
      if(doubles(hI)%i.eq.i) then
       write(170+Nind,'(G20.12,4I5)') doubles(hI)%x,doubles(hI)%i,&
                        doubles(hI)%a,doubles(hI)%j,doubles(hI)%b
        endif
      enddo
    enddo
    close(170+Nind)

    open (unit=170+Nind,file=filecoefTestSrt3,status='old', action='read')
    hI=0
    do
     hI = hI + 1
      read(170+Nind,*,IOSTAT=z) doubles(hI)%x,doubles(hI)%i,&
                  doubles(hI)%a,doubles(hI)%j,doubles(hI)%b
      if (z<0) exit
    enddo
    close(170+Nind)
    do j = 2,nel
     do hI = 1,totex_coeff(2)
      if(doubles(hI)%j.eq.j) then
       write(180+Nind,'(G20.12,4I5)') doubles(hI)%x,doubles(hI)%i,&
                        doubles(hI)%a,doubles(hI)%j,doubles(hI)%b
        endif
      enddo
    enddo
    close(180+Nind)



   ! sorting triples
   else if(Nind.eq.3) then
    do c = 3,nbasis-nel
     do hI = 1,totex_coeff(3)
      if(triples(hI)%c.eq.c.and..not.near_zero(triples(hI)%x)) then
       write(150+Nind,'(G20.12,6I5)') triples(hI)%x,triples(hI)%i,&
                        triples(hI)%a,triples(hI)%j,triples(hI)%b,&
                                      triples(hI)%k,triples(hI)%c
      endif
     enddo
    enddo
    close(150+Nind)

    open (unit=150+Nind,file=filecoefTestSrt1,status='old', action='read')
    hI=0
    do
     hI = hI + 1
     read(150+Nind,*,IOSTAT=z) triples(hI)%x,triples(hI)%i,&
                 triples(hI)%a,triples(hI)%j,triples(hI)%b,&
                               triples(hI)%k,triples(hI)%c
     if (z<0) exit
    enddo
    close(150+Nind)
    do b = 2,nbasis-nel-1
     do hI = 1,totex_coeff(3)
      if(triples(hI)%b.eq.b) then
       write(160+Nind,'(G20.12,6I5)') triples(hI)%x,triples(hI)%i,&
                        triples(hI)%a,triples(hI)%j,triples(hI)%b,&
                                      triples(hI)%k,triples(hI)%c
      endif
     enddo
    enddo
    close(160+Nind)

    open (unit=160+Nind,file=filecoefTestSrt2,status='old', action='read')
    hI=0
    do
     hI = hI + 1
     read(160+Nind,*,IOSTAT=z) triples(hI)%x,triples(hI)%i,&
                 triples(hI)%a,triples(hI)%j,triples(hI)%b,&
                               triples(hI)%k,triples(hI)%c
     if (z<0) exit
    enddo
    close(160+Nind)
    do a = 1,nbasis-nel-2
     do hI = 1,totex_coeff(3)
      if(triples(hI)%a.eq.a) then
       write(170+Nind,'(G20.12,6I5)') triples(hI)%x,triples(hI)%i,&
                        triples(hI)%a,triples(hI)%j,triples(hI)%b,&
                                      triples(hI)%k,triples(hI)%c
      endif
     enddo
    enddo
    close(170+Nind)

    open (unit=170+Nind,file=filecoefTestSrt3,status='old', action='read')
    hI=0
    do
     hI = hI + 1
     read(170+Nind,*,IOSTAT=z) triples(hI)%x,triples(hI)%i,&
                 triples(hI)%a,triples(hI)%j,triples(hI)%b,&
                               triples(hI)%k,triples(hI)%c
     if (z<0) exit
    enddo
    close(170+Nind)
    do i = 1,nel-2
     do hI = 1,totex_coeff(3)
      if(triples(hI)%i.eq.i) then
       write(180+Nind,'(G20.12,6I5)') triples(hI)%x,triples(hI)%i,&
                        triples(hI)%a,triples(hI)%j,triples(hI)%b,&
                                      triples(hI)%k,triples(hI)%c
       endif
      enddo
    enddo
    close(180+Nind)

    open (unit=180+Nind,file=filecoefTestSrt4,status='old', action='read')
    hI=0
    do
     hI = hI + 1
     read(180+Nind,*,IOSTAT=z) triples(hI)%x,triples(hI)%i,&
                 triples(hI)%a,triples(hI)%j,triples(hI)%b,&
                               triples(hI)%k,triples(hI)%c
     if (z<0) exit
    enddo
    close(180+Nind)
    do j = 2,nel-1
     do hI = 1,totex_coeff(3)
      if(triples(hI)%j.eq.j) then
       write(190+Nind,'(G20.12,6I5)') triples(hI)%x,triples(hI)%i,&
                        triples(hI)%a,triples(hI)%j,triples(hI)%b,&
                                      triples(hI)%k,triples(hI)%c
      endif
     enddo
    enddo
    close(190+Nind)

    open (unit=190+Nind,file=filecoefTestSrt5,status='old', action='read')
    hI=0
    do
     hI = hI + 1
     read(190+Nind,*,IOSTAT=z) triples(hI)%x,triples(hI)%i,&
                 triples(hI)%a,triples(hI)%j,triples(hI)%b,&
                               triples(hI)%k,triples(hI)%c
     if (z<0) exit
    enddo
    close(190+Nind)
    do k = 3,nel
     do hI = 1,totex_coeff(3)
      if(triples(hI)%k.eq.k) then
       write(200+Nind,'(G20.12,6I5)') triples(hI)%x,triples(hI)%i,&
                        triples(hI)%a,triples(hI)%j,triples(hI)%b,&
                                      triples(hI)%k,triples(hI)%c
      endif
     enddo
    enddo
    close(200+Nind)
   endif
  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SOLUZIONE IN CUI USO DUE ARRAY, coeffx*(:,:) e coeffx_tmp(:,:),
!! SU CUI SCRIVO E LEGGO VARIE VOLTE FINO A CHE NON È SORTED COME VOGLIO
!!  -> NECI dà problemi di memoria se alloca un 2ndo array for triples
!!
!!    integer :: i, z, Nind, h1, h2, h3, hI
!!    real(dp),allocatable :: coeffx1(:,:),coeffx2(:,:),coeffx3(:,:)
!!    real(dp),allocatable :: coeffx1_tmp(:,:),coeffx2_tmp(:,:),coeffx3_tmp(:,:)
!!    integer :: j,k,a,b,c, ii
!!    integer :: a,b,j 
!!    real(dp):: x,y,l,k,m
!!
!!    allocate(coeffx1(totex_coeff(1),3))
!!    allocate(coeffx2(totex_coeff(2),5))
!!    if(n_store_ci_level.eq.3) allocate(coeffx3(totex_coeff(3),7))
!!
!!    h1=0
!!    h2=0
!!    h3=0
!!    ! loop to read all the CI coefficients from ASCII files
!!    do Nind=1,n_store_ci_level
!!      write(fileCICoeffTest, '( "ci_coeff_",I1,"_test" )') Nind
!!      open (unit=40+Nind,file=fileCICoeffTest,status='old', action='read')
!!      write(filecoefTestSrt0, '("ci_coeff_",I1,"_sort0" )') Nind
!!      open (unit=140+Nind,file=filecoefTestSrt0,status='replace')
!!      write(filecoefTestSrt1, '("ci_coeff_",I1,"_sort1" )') Nind
!!      open (unit=150+Nind,file=filecoefTestSrt1,status='replace')
!!
!!      do
!!        if(Nind.eq.1) then
!!          h1 = h1 + 1
!!          read(40+Nind,*,IOSTAT=z) (coeffx1(h1,i),i=1,3)
!!          if (z<0) then
!!             exit
!!          else if(near_zero(coeffx1(h1,1))) then
!!             cycle
!!          endif
!!        endif
!!        if(Nind.eq.2) then
!!          h2 = h2 + 1
!!          read(40+Nind,*,IOSTAT=z) (coeffx2(h2,i),i=1,5)
!!          if (z<0) then
!!             exit
!!          else if(near_zero(coeffx2(h2,1))) then
!!             cycle
!!          endif
!!        endif
!!        if(Nind.eq.3) then
!!          h3 = h3 + 1
!!          read(40+Nind,*,IOSTAT=z) (coeffx3(h3,i),i=1,7)
!!          if (z<0) then
!!             exit
!!          else if(near_zero(coeffx3(h3,1))) then
!!             cycle
!!          endif
!!        endif
!!      enddo
!!      close (40+Nind)
!!    enddo
!!
!!    ! writing CI coefficients
!!    do hI = 1,totex_coeff(1)
!!       write(141,'(G20.12,2F6.0)') coeffx1(hI,1),coeffx1(hI,2),coeffx1(hI,3)
!!    enddo
!!    close (141)
!!    do hI = 1,totex_coeff(2)
!!       write(142,'(G20.12,4F6.0)') coeffx2(hI,1),coeffx2(hI,2),coeffx2(hI,3),&
!!                                   coeffx2(hI,4),coeffx2(hI,5)
!!    enddo
!!    close (142)
!!    do hI = 1,totex_coeff(3)
!!       write(143,'(G20.12,6F6.0)') coeffx3(hI,1),coeffx3(hI,2),coeffx3(hI,3),&
!!                     coeffx3(hI,4),coeffx3(hI,5),coeffx3(hI,6),coeffx3(hI,7)
!!    enddo
!!    close (143)
!!
!!    ! sorting singles
!!    allocate(coeffx1_tmp(totex_coeff(1),3))
!!    ii=0
!!    do a = 1,nbasis-nel
!!      do h1 = 1,totex_coeff(1)
!!        if(coeffx1(h1,3).eq.a) then
!!          ii = ii + 1
!!          coeffx1_tmp(ii,:) = coeffx1(h2,:)
!!        endif
!!      enddo
!!    enddo
!!    do i = 1,nel
!!      do h1 = 1,totex_coeff(1)
!!        if(coeffx1_tmp(h1,3).eq.i) then
!!          write(151,'(G20.12,2F6.0)') coeffx1_tmp(h1,1),coeffx1_tmp(h1,2),coeffx1_tmp(h1,3)
!!        endif
!!      enddo
!!    enddo
!!    close (151)
!!    deallocate(coeffx1_tmp)
!!
!!    ! sorting doubles
!!    allocate(coeffx2_tmp(totex_coeff(2),5))
!!    ii=0
!!    do a = 1,nbasis-nel-1
!!      do h2 = 1,totex_coeff(2)
!!        if(coeffx2(h2,3).eq.a) then
!!          ii = ii + 1
!!          coeffx2_tmp(ii,:) = coeffx2(h2,:)
!!        endif
!!      enddo
!!    enddo
!!    ii=0
!!    do b = 2,nbasis-nel
!!      do h2 = 1,totex_coeff(2)
!!        if(coeffx2_tmp(h2,5).eq.b) then
!!          ii = ii + 1
!!          coeffx2(ii,:) = coeffx2_tmp(h2,:)
!!        endif
!!      enddo
!!    enddo
!!    ii=0
!!    do i = 1,nel-1
!!      do h2 = 1,totex_coeff(2)
!!        if(coeffx2(h2,2).eq.i) then
!!          ii = ii + 1
!!          coeffx2_tmp(ii,:) = coeffx2(h2,:)
!!        endif
!!      enddo
!!    enddo
!!    do j = 2,nel
!!      do h2 = 1,totex_coeff(2)
!!        if(coeffx2_tmp(h2,4).eq.j) then
!!          write(152,'(G20.12,4F6.0)') coeffx2_tmp(h2,1),coeffx2_tmp(h2,2),coeffx2_tmp(h2,3),&
!!                                      coeffx2_tmp(h2,4),coeffx2_tmp(h2,5)
!!        endif
!!      enddo
!!    enddo
!!    close (152)
!!    deallocate(coeffx2_tmp)
!!
!!    ! sorting triples
!!    allocate(coeffx3_tmp(totex_coeff(3),7))
!!    ii=0
!!    do c = 3,nbasis-nel
!!      do h2 = 1,totex_coeff(3)
!!        if(coeffx3(h2,7).eq.c) then
!!          ii = ii + 1
!!          coeffx3_tmp(ii,:) = coeffx3(h2,:)
!!        endif
!!      enddo
!!    enddo
!!    ii=0
!!    do b = 2,nbasis-nel-1
!!      do h2 = 1,totex_coeff(3)
!!        if(coeffx3_tmp(h2,5).eq.b) then
!!          ii = ii + 1
!!          coeffx3(ii,:) = coeffx3_tmp(h2,:)
!!        endif
!!      enddo
!!    enddo
!!    ii=0
!!    do a = 1,nbasis-nel-2
!!      do h2 = 1,totex_coeff(3)
!!        if(coeffx3(h2,3).eq.a) then
!!          ii = ii + 1
!!          coeffx3_tmp(ii,:) = coeffx3(h2,:)
!!        endif
!!      enddo
!!    enddo
!!    ii=0
!!    do i = 1,nel-2
!!      do h2 = 1,totex_coeff(3)
!!        if(coeffx3(h2,2).eq.i) then
!!          ii = ii + 1
!!          coeffx3(ii,:) = coeffx3_tmp(h2,:)
!!        endif
!!      enddo
!!    enddo
!!    do j = 2,nel-1
!!      do h2 = 1,totex_coeff(3)
!!        if(coeffx3_tmp(h2,4).eq.j) then
!!          ii = ii + 1
!!          coeffx3_tmp(ii,:) = coeffx3(h2,:)
!!        endif
!!      enddo
!!    enddo
!!    do k = 3,nel
!!      do h2 = 1,totex_coeff(3)
!!        if(coeffx3_tmp(h2,6).eq.k) then
!!          write(153,'(G20.12,6F6.0)') coeffx3_tmp(h2,1),coeffx3_tmp(h2,2),&
!!                    coeffx3_tmp(h2,3),coeffx3_tmp(h2,4),coeffx3_tmp(h2,5),&
!!                                      coeffx3_tmp(h2,6),coeffx3_tmp(h2,7)
!!        endif
!!      enddo
!!    enddo
!!    close (153)
!!    deallocate(coeffx3_tmp)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  SOLUZIONE PIU CONVENIENTE SE FUNZIONA
!    ! sorting CI coefficients
!    call sort(coeffx1(:,3))
!    call sort(coeffx1(:,2))
!    call sort(coeffx2(:,3))
!    call sort(coeffx2(:,5))
!    call sort(coeffx2(:,2))
!    call sort(coeffx2(:,4))
!    call sort(coeffx3(:,7))
!    call sort(coeffx3(:,5))
!    call sort(coeffx3(:,3))
!    call sort(coeffx3(:,2))
!    call sort(coeffx3(:,4))
!    call sort(coeffx3(:,6))
!    ! writing sorted CI coefficients
!    do h1 = 1,totex_coeff(1)
!       write(151,'(G20.12,2F6.0)') coeffx1(h1,1),coeffx1(h1,2),coeffx1(h1,3)
!    enddo
!    close (151)
!    do h2 = 1,totex_coeff(2)
!       write(152,'(G20.12,4F6.0)') coeffx2(h2,1),coeffx2(h2,2),coeffx2(h2,3),&
!                                   coeffx2(h2,4),coeffx2(h2,5)
!    enddo
!    close (152)
!    do h3 = 1,totex_coeff(3)
!       write(153,'(G20.12,6F6.0)') coeffx3(h3,1),coeffx3(h3,2),coeffx3(h3,3),&
!                     coeffx3(h3,4),coeffx3(h3,5),coeffx3(h3,6),coeffx3(h3,7)
!    enddo
!    close (153)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! OLD SORTING
!   ! Writing the CI coefficients in the ASCII files
!    do i = 1,nel
!      do a = nel+1,nbasis
!        if(.not. near_zero(S(i,a))) then
!          write(101,'(G20.12,2I5)') S(i,a),i,a-nel
!        endif
!      enddo
!    enddo
!    close(101)
!
!    do j = 2,nel
!      do i = 1, nel-1
!        do b = nel+1, nbasis
!          do a = nel,nbasis-1
!            if(.not. near_zero(D(i,a,j,b))) then
!              write(102,'(G20.12,4I5)') D(i,a,j,b),i,a-nel,j,b-nel
!            endif
!          enddo
!        enddo
!      enddo
!    enddo
!    close(102)
!
!    if(n_store_ci_level.eq.3) then
!      do k = 3,nel
!        do j = 2, nel-1
!          do i = 1, nel-2
!            do a = nel+1,nbasis
!              do b = a+1, nbasis
!                do c = b+1, nbasis
!                  if(.not. near_zero(T(i,a,j,b,k,c))) then
!                    write(103,'(G20.12,6I5)') T(i,a,j,b,k,c),i,a-nel,j,b-nel,k,c-nel
!                  endif
!                enddo
!              enddo
!            enddo
!          enddo
!        enddo
!      enddo
!      close(103)
!    endif

  end subroutine sorting

end module sdt_amplitudes
