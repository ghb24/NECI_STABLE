#:include "macros.fpph"
#:include "algorithms.fpph"

! Module to collect and average the CI coefficients
module sdt_amplitudes

  use bit_reps, only: extract_sign, decode_bit_det, encode_sign, niftot, nifd
  use constants, only: dp, lenof_sign, n_int, int64, stdout
  use DetBitOps, only: get_bit_excitmat, EncodeBitDet, GetBitExcitation!, FindBitExcitLevel
  use util_mod, only: near_zero
  use FciMCData, only: TotWalkers, iLutRef, CurrentDets, AllNoatHF, projedet, &
                       ll_node
  use hash, only: hash_table_lookup, add_hash_table_entry, init_hash_table, &
                  clear_hash_table
  use LoggingData, only: n_store_ci_level, n_iter_ci_coeff
  use SystemData, only: nel, nbasis, symmax
  use Parallel_neci, only: iProcIndex, MPIcollection
  use MPI_wrapper, only: root
  use sort_mod, only: sort

  implicit none

  type :: t_singles
    double precision :: x
    integer :: i,a
  end type
  type :: t_doubles
    double precision :: x
    integer :: i,a,j,b
  end type
  type :: t_triples
    double precision :: x
    integer :: i,a,j,b,k,c
  end type

  integer(n_int), allocatable :: ciCoeff_storage(:,:), root_ciCoeff_storage(:,:)
  integer :: hash_table_ciCoeff_size, first_free_entry, nCyc, root_first_free_entry
  type(ll_node), pointer :: hash_table_ciCoeff(:)
  character (len=90) :: fileCICoeffAv, fileCIcoeffSort
  integer(n_int), allocatable  :: totEntCoeff(:,:)

contains

   logical elemental function sing_a(p1, p2)
       type(t_singles), intent(in) :: p1, p2
       sing_a = p1%a <= p2%a
   end function
   logical elemental function sing_i(p1, p2)
       type(t_singles), intent(in) :: p1, p2
       sing_i = p1%i <= p2%i
   end function

   logical elemental function doub_a(p1, p2)
       type(t_doubles), intent(in) :: p1, p2
       doub_a = p1%a <= p2%a
   end function
   logical elemental function doub_b(p1, p2)
       type(t_doubles), intent(in) :: p1, p2
       doub_b = p1%b <= p2%b
   end function
   logical elemental function doub_i(p1, p2)
       type(t_doubles), intent(in) :: p1, p2
       doub_i = p1%i <= p2%i
   end function
   logical elemental function doub_j(p1, p2)
       type(t_doubles), intent(in) :: p1, p2
       doub_j = p1%j <= p2%j
   end function

   logical elemental function trip_a(p1, p2)
       type(t_triples), intent(in) :: p1, p2
       trip_a = p1%a <= p2%a
   end function
   logical elemental function trip_b(p1, p2)
       type(t_triples), intent(in) :: p1, p2
       trip_b = p1%b <= p2%b
   end function
   logical elemental function trip_c(p1, p2)
       type(t_triples), intent(in) :: p1, p2
       trip_c = p1%c <= p2%c
   end function
   logical elemental function trip_i(p1, p2)
       type(t_triples), intent(in) :: p1, p2
       trip_i = p1%i <= p2%i
   end function
   logical elemental function trip_j(p1, p2)
       type(t_triples), intent(in) :: p1, p2
       trip_j = p1%j <= p2%j
   end function
   logical elemental function trip_k(p1, p2)
       type(t_triples), intent(in) :: p1, p2
       trip_k = p1%k <= p2%k
   end function

  subroutine init_ciCoeff
    nCyc = 0
    first_free_entry = 0
    hash_table_ciCoeff_size = 500000
    allocate(hash_table_ciCoeff(hash_table_ciCoeff_size))
    call init_hash_table(hash_table_ciCoeff)
    allocate(ciCoeff_storage(0:NIfTot,hash_table_ciCoeff_size))
    allocate(totEntCoeff(n_store_ci_level,2))
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

    integer  :: i,ic,ex(2,4),icI
    real(dp) :: sign_tmp(lenof_sign), ref_coef
    logical  :: tPar

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

      totEntCoeff(:,:) = 0
      do icI = 0, n_store_ci_level
         if(icI.ne.0) then
           write(fileCICoeffAv, '( "ci_coeff_",I1,"_av" )') icI
           open (unit=30+icI,file=fileCICoeffAv,status='replace')
         endif
         ! loop over the total entries of CI coefficients
         do i = 1, root_first_free_entry
            ic = n_store_ci_level+1
            ! gets the excitation level of the CI coefficient
            call get_bit_excitmat(iLutRef(:,1),root_ciCoeff_storage(:,i),ex,ic)
!            ic = FindBitExcitLevel(ilutRef(:,1),root_ciCoeff_storage(:,i))
            if(icI.eq.ic) then
              ! gets the value of the CI coefficient (i.e. the number of walkers)
              call extract_sign(root_ciCoeff_storage(:,i),sign_tmp)
              ex(1,1) = ic
              ! gets the sign of the CI coef (tPar=true if odd number of permutations)
              call GetBitExcitation(ilutRef(:,1),root_ciCoeff_storage(:,i),ex,tPar)
              if(tPar) sign_tmp = -sign_tmp

               select case(icI) ! writing averaged CI coefficients
               case(0)  ! reference
                  write(stdout,"(A44,F14.3)") 'Instantaneous number of walkers on HF     = ', AllNoatHF
                  ref_coef = -sign_tmp(1)
                  write(stdout,"(A44,F14.3)") 'Averaged number of walkers on HF          = ', -sign_tmp/nCyc
                  write(stdout,"(A44,I10)") 'Total entries of CI coefficients          = ', root_first_free_entry
               case(1)  ! singles
                 totEntCoeff(icI,1) = totEntCoeff(icI,1) + 1        ! total entries for singles
                 if(.not. near_zero(sign_tmp(1))) then
                   totEntCoeff(icI,2) = totEntCoeff(icI,2) + 1      ! total entries for singles without zeros
                   write(31,'(G20.12,2I5)') sign_tmp/ref_coef,ex(1,1),ex(2,1)
                 endif
               case(2)  ! doubles
                 totEntCoeff(icI,1) = totEntCoeff(icI,1) + 1        ! total entries for doubles
                 if(.not. near_zero(sign_tmp(1))) then
                   totEntCoeff(icI,2) = totEntCoeff(icI,2) + 1      ! total entries for doubles without zeros
                   write(32,'(G20.12,4I5)') sign_tmp/ref_coef,ex(1,1),&
                                            ex(2,1),ex(1,2),ex(2,2)
                 endif
               case(3)  ! triples
                 totEntCoeff(icI,1) = totEntCoeff(icI,1) + 1        ! total entries for triples
                 if(.not. near_zero(sign_tmp(1))) then
                   totEntCoeff(icI,2) = totEntCoeff(icI,2) + 1      ! total entries for triples without zeros
                   write(33,'(G20.12,6I5)') sign_tmp/ref_coef,ex(1,1),&
                                            ex(2,1),ex(1,2),ex(2,2),ex(1,3),ex(2,3)
                 endif
               end select
            end if
         enddo
         if(iCI.ne.0) then
           close(30+icI)
           write(stdout,"(A28,I1,A14,I11)") ' total entries for ci_coeff_',icI,' =', totEntCoeff(iCI,1)
           if(totEntCoeff(iCI,1).ne.totEntCoeff(iCI,2)) then
             write(stdout,"(A25,I1,A17,I11)") ' -total entries ci_coeff_',icI,&
                                              ' without zeros  =', totEntCoeff(iCI,2)
           endif
         endif
      enddo

      write(stdout,*) ' sorting CI coefficients...'
      call sorting()

      write(stdout,*) '-> CI coefficients written in ASCII files ci_coeff_*'
      write(stdout,*) ''
      write(stdout,*) '-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --'
    endif

    call fin_ciCoeff()
  end subroutine print_averaged_ci_coeff


  ! it lists the averaged CI coeffs sorting the indices in
  ! this way: OCC(alpha), OCC(beta), VIR(alpha), VIR(beta)
 subroutine sorting

  integer  :: i,hI,iCI,signCI,ex(2,n_store_ci_level),Itot(nbasis)
  real(dp) :: x
  type(t_singles),allocatable :: singles(:)
  type(t_doubles),allocatable :: doubles(:)
  type(t_triples),allocatable :: triples(:)


  call alphaBetaOrbs(Itot)

  do iCI=1,n_store_ci_level
    write(fileCICoeffAv, '( "ci_coeff_",I1,"_av" )') iCI
    open (unit=30+iCI,file=fileCICoeffAv,status='old', action='read')
    write(fileCIcoeffSort, '("ci_coeff_",I1 )') iCI
    open (unit=140+iCI,file=fileCIcoeffSort,status='replace')

    if(iCI.eq.1) allocate(singles(totEntCoeff(iCI,2)))
    if(iCI.eq.2) allocate(doubles(totEntCoeff(iCI,2)))
    if(iCI.eq.3) allocate(triples(totEntCoeff(iCI,2)))

    do hI=1,totEntCoeff(iCI,2)
      read(30+iCI,*) x,(ex(1,i),ex(2,i),i=1,icI)

      call indPreSort(icI,Itot,ex,signCI)

      if(iCI.eq.1) then
        singles(hI)%x = signCI*x
        singles(hI)%i = ex(1,1)
        singles(hI)%a = ex(2,1)
      else if(iCI.eq.2) then
        doubles(hI)%x = signCI*x
        doubles(hI)%i = ex(1,1)
        doubles(hI)%a = ex(2,1)
        doubles(hI)%j = ex(1,2)
        doubles(hI)%b = ex(2,2)
      else if(iCI.eq.3) then
        triples(hI)%x = signCI*x
        triples(hI)%i = ex(1,1)
        triples(hI)%a = ex(2,1)
        triples(hI)%j = ex(1,2)
        triples(hI)%b = ex(2,2)
        triples(hI)%k = ex(1,3)
        triples(hI)%c = ex(2,3)
      endif
    enddo

    select case(iCI) ! sorting and writing CI coefficients
    case(1)  ! singles
      @:sort(t_singles, singles, rank=1, along=1, comp=sing_a)
      @:sort(t_singles, singles, rank=1, along=1, comp=sing_i)
      do hI = 1,totEntCoeff(iCI,2)
        write(140+iCI,'(G20.12,2I5)') singles(hI)%x,singles(hI)%i,&
                                      singles(hI)%a
      enddo
      close (140+iCI)
      deallocate(singles)

    case(2)  ! doubles
      @:sort(t_doubles, doubles, rank=1, along=1, comp=doub_a)
      @:sort(t_doubles, doubles, rank=1, along=1, comp=doub_b)
      @:sort(t_doubles, doubles, rank=1, along=1, comp=doub_i)
      @:sort(t_doubles, doubles, rank=1, along=1, comp=doub_j)
      do hI = 1,totEntCoeff(iCI,2)
         write(140+iCI,'(G20.12,4I5)') doubles(hI)%x,doubles(hI)%i,&
                         doubles(hI)%a,doubles(hI)%j,doubles(hI)%b
      enddo
      close (140+iCI)
      deallocate(doubles)

    case(3)  ! triples
      @:sort(t_triples, triples, rank=1, along=1, comp=trip_c)
      @:sort(t_triples, triples, rank=1, along=1, comp=trip_b)
      @:sort(t_triples, triples, rank=1, along=1, comp=trip_a)
      @:sort(t_triples, triples, rank=1, along=1, comp=trip_i)
      @:sort(t_triples, triples, rank=1, along=1, comp=trip_j)
      @:sort(t_triples, triples, rank=1, along=1, comp=trip_k)
      do hI = 1,totEntCoeff(iCI,2)
        write(140+iCI,'(G20.12,6I5)') triples(hI)%x,triples(hI)%i,&
                        triples(hI)%a,triples(hI)%j,triples(hI)%b,&
                                      triples(hI)%k,triples(hI)%c
      enddo
      close (140+iCI)
      deallocate(triples)
    end select

  enddo ! do iCI=1,n_store_ci_level

 end subroutine sorting


  ! it finds all the alpha/beta occ/unocc orbs
  subroutine alphaBetaOrbs(Itot)

    use SymExcitDataMod, only:  OrbClassCount
    use GenRandSymExcitNUMod , only : ClassCountInd

    integer, intent(out) :: Itot(nbasis)
    integer  :: i,j,k,z,ial(symmax),ibe(symmax),iSym,totEl,totOrb
    integer  :: Norb(symmax),NorbTot(symmax)

      NorbTot(:)=0
      i=1
      do iSym = 1, symmax
        i = i + 1
        Norb(iSym) =  OrbClassCount(ClassCountInd(1/2,iSym,0)) * 2
        NorbTot(i) = Norb(iSym) + NorbTot(i-1)
      enddo
     ! loop to find all the occupied alpha spin orbitals
      totEl = 0
      do iSym = 1, symmax
        ial(iSym)=0
        do z=1,nel
          if(projEDet(z,1).gt.NorbTot(iSym) .and. projEDet(z,1).le.NorbTot(iSym+1)) then
            if(MOD(projEDet(z,1),2).eq.1) then
              ial(iSym)=ial(iSym)+1
              Itot(totEl+ial(iSym))=projEDet(z,1)
            endif
          endif
        enddo
        totEl = totEl + ial(iSym)
        ibe(iSym)=0
     ! loop to find all the occupied beta spin orbitals
        do z=1,nel
          if(projEDet(z,1).gt.NorbTot(iSym) .and. projEDet(z,1).le.NorbTot(iSym+1)) then
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
        do k=NorbTot(iSym)+1,NorbTot(iSym+1)
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
        do k=NorbTot(iSym)+1,NorbTot(iSym+1)
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

  end subroutine alphaBetaOrbs


  ! PreSorting routine which converts the indices into Molpro standard
  subroutine indPreSort(icI,Itot,ex,signCI)

    use util_mod, only: swap

    integer, intent(in)  :: icI,Itot(nbasis)
    integer, intent(out) :: signCI
    integer, intent(inout) :: ex(2,n_store_ci_level)
    integer  :: j,k,p

    ! indices conversion
    do k=1,icI
      do p=1,nel
          if(ex(1,k).eq.Itot(p)) then
            ex(1,k)=p
            exit
         endif
      enddo
      do p=nel+1,nbasis
         if(ex(2,k).eq.Itot(p)) then
            ex(2,k)=p-nel
            exit
          endif
      enddo
    enddo
    ! insertion sort
    signCI=1
    do p=1,2
      do k=2,icI
        do j=k,2,-1
          if (ex(p,j).lt.ex(p,j-1)) then
            call swap(ex(p,j),ex(p,j-1))
            signCI=-signCI
          else
            exit
          endif
        enddo
      enddo
    enddo

  end subroutine indPreSort


  subroutine fin_ciCoeff
    call clear_hash_table(hash_table_ciCoeff)
    deallocate(hash_table_ciCoeff)
    deallocate(ciCoeff_storage)
    deallocate(totEntCoeff)
  end subroutine fin_ciCoeff

end module sdt_amplitudes
