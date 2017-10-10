#include "macros.h"
module adi_references
use Parallel_neci
use FciMCData, only: ilutRefAdi, ilutRef, nRefs, nIRef, signsRef, nRefsSings, nRefsDoubs
use bit_rep_data, only: niftot, nifdbo
use constants
use SystemData, only: nel

implicit none

contains

  subroutine setup_reference_space(tPopPresent)
    use CalcData, only: tReadRefs, tProductReferences, nExProd
    use LoggingData, only: ref_filename
    implicit none
    logical, intent(in) :: tPopPresent
    integer :: nRead, nRCOld
    logical :: tGen

    ! A first guess at the number of references
    nRefs = max(nRefsDoubs, nRefsSings)

    call reallocate_ilutRefAdi(nRefs)

    ! If we actually did something
    tGen = .false.
    ! First, generate the reference space of nRefs determinants from the population (if present)
    if(tPopPresent) then
       call generate_ref_space()
       tGen = .true.
    endif

    if(tReadRefs) then
       ! Then, add the references from the file
       call read_in_refs(ref_filename, nRead, tPopPresent)
       ! If we added references, note this
       if(nRead > nRefs) nRefs = nRead
       tGen = .true.
    endif
    ! Then, add the product references
    if(tGen) then
       nRCOld = nRefs
       if(tProductReferences) then
          call add_product_references(nRefs, nExProd)
          call update_ref_signs()
       endif
       ! And prompt the output message
       call print_reference_notification(nRefs, nRCOld)
    endif
    
  end subroutine setup_reference_space

!------------------------------------------------------------------------------------------!

  subroutine read_in_refs(filename, nRead, tPopPresent)
    use util_mod, only: get_free_unit
    use CalcData, only: tDelayGetRefs
    implicit none
    integer, intent(out) :: nRead
    character(255), intent(in) :: filename
    logical, intent(in) :: tPopPresent
    integer :: iunit, i, run, stat
    integer(n_int) :: tmp(0:NIfTot)
    logical :: exists
    character(*), parameter :: this_routine = "read_in_refs"

    ! All procs read in the file, as all need all references
    inquire(file=trim(filename), exist = exists)
    if(.not. exists) call stop_all(this_routine, "No "//trim(filename)//" file detected.")
    iunit = get_free_unit()

    ! If nRefs == 1, we set it to the number of references in the file
    if(nRefs == 1) then
       ! Therefore, we scan the file once
       open(iunit, file=trim(filename), status = 'old')
       nRefs = 0
       do
          read(iunit, *, iostat = stat) tmp
          if(stat < 0) exit
          nRefs = nRefs + 1
       end do
       close(iunit)

       ! And allocate the ilutRefAdi accordingly
       call reallocate_ilutRefAdi(nRefs)
    endif

    open(iunit, file=trim(filename), status = 'old')

    ! Read in at most nRefs references
    nRead = 0
    do i = 1, nRefs
       ! Note that read-in references always have precedence over generated references
       read(iunit, *, iostat=stat) ilutRefAdi(:,1,i)
       ! If there are no more dets to be read, exit
       if(stat < 0) exit
       ! If we successfully read, log it
       nRead = nRead + 1
    enddo

    ! Copy the references to all runs
    do i = 1, nRefs
       do run = 2, inum_runs
          ilutRefAdi(:,run,i) = ilutRefAdi(:,1,i)
       end do
    enddo

    close(iunit)
    ! If we read in less than nRefs detererminants, we want to get more later
    if(nRead < max(nRefsSings, nRefsDoubs) .and. .not. tPopPresent) tDelayGetRefs = .true.
  end subroutine read_in_refs

!------------------------------------------------------------------------------------------!

    subroutine generate_ref_space()
      use semi_stoch_gen, only: generate_space_most_populated
      use LoggingData, only: ref_filename, tWriteRefs
      implicit none
      integer(MPIArg) :: mpi_refs_found
      integer :: ierr, i, all_refs_found, refs_found
      integer(MPIArg) :: refs_found_per_proc(0:nProcessors-1), refs_displs(0:nProcessors-1)
      integer(n_int) :: ref_buf(0:NIfTot,nRefs), mpi_buf(0:NIfTot,nRefs)
      character(*), parameter :: this_routine = "generate_ref_space"
      
      write(6,*) "Getting reference determinants for all-doubs-initiators"

      ! we need to be sure ilutRefAdi has the right size
      call reallocate_ilutRefAdi(nRefs)

      ! Get the nRefs most populated determinants
      refs_found = 0
      call generate_space_most_populated(nRefs, .false., 1, ref_buf, refs_found)
      ! Communicate the refs_found info
      mpi_refs_found = int(refs_found,MPIArg)
      call MPIAllGather(mpi_refs_found, refs_found_per_proc, ierr)
      all_refs_found = sum(refs_found_per_proc)
      if(all_refs_found .ne. nRefs) then
         write(6,*) "all_refs_found = ", all_refs_found
         call stop_all(this_routine, &
           "Number of references found differs from requested")
      endif
      refs_displs(0) = 0
      do i = 1, nProcessors - 1
         refs_displs(i) = sum(refs_found_per_proc(0:i-1))
      enddo
      ! Store them on all processors
      call MPIAllGatherV(ref_buf(0:NIfTot, 1:refs_found), mpi_buf, refs_found_per_proc, refs_displs)

      ! And write the so-merged references to ilutRef
      call set_additional_references(mpi_buf)

      if(tWriteRefs) call output_reference_space(ref_filename)
    end subroutine generate_ref_space

!------------------------------------------------------------------------------------------!

    subroutine set_additional_references(mpi_buf)
      use DetBitOps, only: DetBitEQ, sign_lt, sign_gt
      use sort_mod, only: sort
      use bit_rep_data, only: extract_sign
      implicit none
      integer(n_int), intent(inout) :: mpi_buf(0:NIfTot,nRefs)
      integer :: i, k, run
      integer(n_int) :: tmp_ilut(0:NIfTot)
      real(dp) :: minOcc, tmp_sgn(lenof_sign)
      logical :: get_sign

      get_sign = .false.
      ! Check if the current reference is in mpi_buf
      ! k is the index of the current reference (-1 if not present)
      k = -1
      do i = 1, nRefs
         if(DetBitEQ(ilutRef(:,1),mpi_buf(:,i),NIfDBO)) k = i
      enddo
      
      ! If the reference is not at pos. 1, move it there
      if(k .ne. -1) then 
         ! If k==1, everything is already alright, nothing to do here
         if(k > 1) then
            tmp_ilut = mpi_buf(:,k)
            mpi_buf(:,k) = mpi_buf(:,1)
            mpi_buf(:,1) = tmp_ilut
         endif
      else
         ! if it is not even in mpi_buf, replace the lowest populated det
         ! in the buffer with it
         k = 1
         call extract_sign(mpi_buf(:,1), tmp_sgn)
         minOcc = sum(abs(tmp_sgn))
         do i = 2, nRefs
            call extract_sign(mpi_buf(:,i),tmp_sgn)
            if(sum(abs(tmp_sgn)) < minOcc) k = i
         enddo
         tmp_ilut = mpi_buf(:,1)
         mpi_buf(:,1) = ilutRef(:,1)
         if(k > 1) mpi_buf(:,k) = tmp_ilut
         ! also, we now have to obtain the sign for the ilutref anew as ilutref does not
         ! contain a sign
         get_sign = .true.
      endif
      ! Now we have the main reference at position 1

      ! Then sort the buffer with respect to sign (but keep the original reference at pos. 1
      call sort(mpi_buf(:,2:), sign_gt, sign_lt)

      ! Copy the buffer to ilutRef
      do i = 1, nRefs
         do run = 1, inum_runs
            ilutRefAdi(:,run,i) = mpi_buf(:,i)
         enddo
      enddo

      ! if needed, update the signs of ilutRefAdi
      if(get_sign) call update_ref_signs()
    end subroutine set_additional_references

!------------------------------------------------------------------------------------------!

    subroutine output_reference_space(filename)
      use ParallelHelper, only: root
      use util_mod, only: get_free_unit
      use FciMCData, only: nRefs
      implicit none
      character(255), intent(in) :: filename
      integer :: iunit, i, ierr
      
      if(iProcIndex == root) then
         iunit = get_free_unit()
         open(iunit, file = filename, status = 'replace')
         ! write the references into the file
         do i = 1, nRefs
            write(iunit, *) ilutRefAdi(:,1,i)
         enddo
         call neci_flush(iunit)
         close(iunit)
      endif
      call MPIBarrier(ierr)
    end subroutine output_reference_space

!------------------------------------------------------------------------------------------!

    subroutine print_reference_notification(nRefs, nRefsPrev)
      use bit_rep_data, only: extract_sign
      implicit none
      integer, intent(in) :: nRefs, nRefsPrev
      integer :: i, j
      real(dp) :: temp_sgn(lenof_sign)
      
      if(iProcIndex==root) then
         print *, "References for all-doubs-initiators set as"
         do i=1, nRefs
            if(i .eq. (nRefsPrev+1)) write(iout,*) "References generated from product excitations:"
            ! First output the reference determinant
            call WriteDetBit(iout,ilutRefAdi(:,1,i),.false.)
            ! And then the sign
            call extract_sign(ilutRefAdi(:,1,i),temp_sgn)
            do j = 1, lenof_sign
               write(iout,"(G16.7)",advance='no') temp_sgn(j)
            enddo
            write(iout,'()')
         enddo
      endif
    end subroutine print_reference_notification

!------------------------------------------------------------------------------------------!

    subroutine update_ref_signs()
      implicit none
      integer :: iRef
      
      do iRef = 1, nRefs
         call update_single_ref_sign(iRef)
      enddo
    end subroutine update_ref_signs

!------------------------------------------------------------------------------------------!

    subroutine update_single_ref_sign(iRef)
      ! Get the signs of the references from the currentdets. Used both in the final
      ! output and in generating the product excitations
      use bit_reps, only: decode_bit_det, encode_sign
      use bit_rep_data, only: NIfDBO, extract_sign
      use hash, only: hash_table_lookup
      use FciMCData, only: CurrentDets, HashIndex
      use Parallel_neci, only: MPISumAll
      implicit none
      integer, intent(in) :: iRef
      integer:: nI(nel), hash_val, index, run
      real(dp) :: tmp_sgn(lenof_sign), mpi_sgn(lenof_sign)
      logical :: tSuccess      
      

      call decode_bit_det(nI,ilutRefAdi(:,1,iRef))
      call hash_table_lookup(nI, ilutRefAdi(:,1,iRef),NIfDBO,HashIndex,CurrentDets,&
           index,hash_val,tSuccess)
      if(tSuccess) then
         call extract_sign(CurrentDets(:,index),tmp_sgn)
      else
         tmp_sgn = 0.0_dp
      endif
      ! This does the trick of communication: We need to get the sign to all 
      ! processors, so just sum up the individual ones
      call MPISumAll(tmp_sgn, mpi_sgn)

      ! The sign contains information on all runs, so it is the same on all
      do run = 1, inum_runs
         call encode_sign(ilutRefAdi(:,run,iRef),mpi_sgn)
      enddo

    end subroutine update_single_ref_sign

!------------------------------------------------------------------------------------------!

    subroutine spin_symmetrize_references()
      use DetBitOps, only: spin_flip, DetBitEQ
      ! If using HPHF, we want to have the spinflipped version of each reference, too
      ! This is a bit cumbersome to implement, as the spinflipped version might or
      ! might not be already in the references
      implicit none
      integer :: iRef, jRef,  cRef, run
      integer(n_int) :: ilutRef_new(0:NIfTot,2*nRefs), tmp_ilut(0:NIfTot)
      logical :: missing
      
      cRef = 1
      do iRef = 1 ,nRefs
         ! First, copy the current ilut to a new buffer
         ilutRef_new(:,cRef) = ilutRefAdi(:,1,iRef)
         cRef = cRef + 1
         ! get the spinflipped determinant
         tmp_ilut = spin_flip(ilutRefAdi(:,1,iRef))
         missing = .true.
         ! check if it is already present
         do jRef = 1, nRefs
            if(DetBitEQ(tmp_ilut,ilutRefAdi(:,1,jRef),NIfDBO)) missing = .false.
         end do
         ! If not, also add it to the list
         if(missing) then
            ilutRef_new(:,cRef) = tmp_ilut
            cRef = cRef + 1
         endif
      enddo

      ! Resize the ilutRef array
      call reallocate_ilutRefAdi(cRef)
      nRefs = cRef
      ! Now, copy the newly constructed references back to the original array
      do iRef = 1, cRef
         do run = 1, inum_runs
            ilutRefAdi(:,run,iRef) = ilutRef_new(:,iRef)
         end do
      end do
    end subroutine spin_symmetrize_references

!------------------------------------------------------------------------------------------!

    subroutine add_product_references(nRefsIn, prodLvl)
      implicit none
      integer, intent(in) :: nRefsIn
      integer, intent(in) :: prodLvl
      integer(n_int) :: ilut_tmp(0:NIfTot)
      integer(n_int), allocatable :: prod_buffer(:,:)
      integer :: nProdsMax, nProds, i, run, cLvl
      logical :: t_is_valid

      ! Get the maximum number of product excitations
      nProdsMax = maxProdEx(nRefsIn, prodLvl)

      ! Setup the buffer
      allocate(prod_buffer(0:NifTot,nProdsMax))
      prod_buffer(:,1:nRefsIn) = ilutRefAdi(:,1,1:nRefsIn)
      nProds = nRefsIn

      do cLvl = 2, prodLvl
         ! For each product level, get all possible product excitations
         do i = 1, nRefsIn**cLvl
            ! Get the excitation product corresponding to i
            call get_product_excitation(i, cLvl, nRefsIn, ilut_tmp, t_is_valid)
            
            ! Check, if it is already in the list (e.g. if both tau operators are those of
            ! the reference)
            if(t_is_valid) t_is_valid = ilut_not_in_list(ilut_tmp, prod_buffer, nProds)
         
            ! If all excitations were valid and the new ilut is not already in the buffer, 
            ! add it
            if(t_is_valid) then
               nProds = nProds + 1
               prod_buffer(:,nProds) = ilut_tmp
            endif
         end do
         ! Do this for each targeted product level
      end do
            
      ! Write the buffer to ilutRefAdi
      call reallocate_ilutRefAdi(nProds)
      ! Reassign the number of references (do not change the actual nRefs, as this
      ! might be important later on)
      nRefs = nProds
      do run = 1, inum_runs
         ilutRefAdi(:,run,:) = prod_buffer(:,1:nProds)
      end do

      deallocate(prod_buffer)
    end subroutine add_product_references

!------------------------------------------------------------------------------------------!

    subroutine get_product_excitation(i, cLvl, nRefs, ilut_tmp, isValid)
      use DetBitOps, only: CountBits
      use bit_rep_data, only: NIfD
      implicit none
      integer, intent(in) :: i, cLvl, nRefs
      integer(n_int), intent(out) :: ilut_tmp(0:NIfTot)
      logical, intent(out) :: isValid
      integer(n_int) :: tau(0:NIfTot), tau_cc(0:NIfTot)
      integer :: cp

      ilut_tmp = ilutRef(:,1)
      isValid = .true.
      tau_cc = 0_n_int
      ! By getting the cLvl excitation operators that correspond to i
      do cp = 1, cLvl
         ! Store the excitation operators in tau
         tau = IEOR(ilutRefAdi(:,1,cpIndex(i,nRefs,cp)),ilutRef(:,1))

         ! Check if we do not annihilate something twice
         if(any(IAND(tau_cc(0:NIfD),tau(0:NIfD)) .ne. 0_n_int)) then
            isValid = .false.
            exit
         endif
         ! And apply it on the temporary iLut
         tau_cc(0:NIfD) = IEOR(tau_cc(0:NIfD),tau(0:NIfD))

      end do

      ! Check if it's valid, if not, go to the next value of i
      if(CountBits(ilut_tmp,NIfD) .ne. nEl) isValid = .false.
      ilut_tmp(0:NIfD) = IEOR(ilut_tmp(0:NIfD), tau_cc(0:NIfD))
    end subroutine get_product_excitation

!------------------------------------------------------------------------------------------!

    subroutine generate_type_one_refs()
      use FciMCData, only: CurrentDets, TotWalkers
      use bit_reps, only: extract_sign
      implicit none
      integer :: i, nTOne, nBlocks, ierr
      integer(n_int), allocatable :: refBuf(:,:)
      real(dp) :: tmp_sign(lenof_sign)
      logical :: is_tone
      integer, parameter :: allocBlock = 100
      character(*), parameter :: this_routine = "generate_type_one_refs"

      allocate(refBuf(0:NIfTot, allocBlock))
      nBlocks = 1
      
      nTOne = 0
      do i = 1, TotWalkers
         call extract_sign(CurrentDets(:,i), tmp_sign)
         ! Check if a determinant should be a type-1 reference
         is_tone = check_type_one_ref(CurrentDets(:,i), tmp_sign)
         if(is_tone) then
            nTOne = nTOne + 1
            ! If we exceed the memory of the buffer, add another block
            if(nTOne > nBlocks*allocBlock) then
               nBlocks = nBlocks + 1
               call resize_ilut_list(refBuf, (nBlocks-1) * allocBlock ,nBlocks * allocBlock)
               if(ierr) call stop_all(this_routine, "Unable to allocate temporary list")
            endif
            ! Add the determinant to the temporary list
            refBuf(:,nTOne) = CurrentDets(:,i)
         endif
      end do

      ! And construct the array used later
      call fill_ilutRefTOne(refBuf, nTOne)
      
    end subroutine generate_type_one_refs

!------------------------------------------------------------------------------------------!

    subroutine fill_ilutRefTOne(list, listSize)
      implicit none
      integer, intent(in) :: listSize
      integer(n_int), intent(in) :: list(0:NIfTot,listSize)
      
      ! resize the ilutRefAdi array
      ! call resize_ilut_list(ilutRefAdi, nRefs, nRefs + listSize)
      ! and add the new iluts
      ilutRefAdi(:,1,(nRefs + 1):(nRefs + listSize)) = list(:,:)
      nRefs = nRefs + listSize
    end subroutine fill_ilutRefTOne


!------------------------------------------------------------------------------------------!

    function check_type_one_ref(ilut, ilut_sign) result(is_tone)
      use DetBitOps, only: FindBitExcitLevel
      implicit none
      integer(n_int), intent(in) :: ilut(0:NIfTot)
      ! We pass the sign as an extra argument to prevent redundant extraction
      real(dp), intent(in) :: ilut_sign(lenof_sign)
      integer :: iRef, exLevel
      logical :: is_coherent
      logical :: is_tone
      
      is_tone = .false.
      is_coherent = .false.
      ! Go through all t-0 references
      do iRef = 1, nRefs
         ! Get the excitation level, only proceed if there can be some coupling
         exLevel = FindBitExcitLevel(ilutRefAdi(:,1,iRef), ilut)
         if(exLevel < 3) then
            ! If we are coherent with two or more, set the t-1 flag
            if(is_coherent) is_tone = .true.
            is_coherent = check_sign_coherence(ilut, ilut_sign, iRef)
            ! If an incoherent t-0 ref is found, return .false.
            if(.not. is_coherent) then
               is_tone = .false.
               return
            endif
         endif
      enddo
    end function check_type_one_ref

!------------------------------------------------------------------------------------------!

    subroutine resize_ilut_list(list, oldSize, newSize)
      implicit none
      integer, intent(in) :: newSize, oldSize
      integer(n_int), allocatable, intent(inout) :: list(:,:)
      integer(n_int) :: tmp(0:NIfTot,oldSize)
      integer :: ierr
      character(*), parameter :: this_routine = "resize_ilut_list"

      ! If the list is already allocated, clear it
      if(allocated(list)) then
         tmp(:,:) = list(:,:)
         deallocate(list)
      endif
      allocate(list(0:NIfTot,newSize), stat = ierr)
      if(ierr) call stop_all(this_routine, "Not enough memory to resize ilut list")
      
      ! Write the old entries of list into the new memory
      list(:,1:oldSize) = tmp(:,1:oldSize)
      
    end subroutine resize_ilut_list

!------------------------------------------------------------------------------------------!

    subroutine reallocate_ilutRefAdi(size)
      implicit none
      integer, intent(in) :: size
      
      ! reallocate first the ilutref itself
      if(allocated(ilutRefAdi)) deallocate(ilutRefAdi)
      allocate(ilutRefAdi(0:NIfTot,inum_runs,size))

      ! and then the corresponding chaches for signs and determinants
      if(allocated(nIRef)) deallocate(nIRef)
      allocate(nIRef(nel,size))
      
      if(allocated(signsRef)) deallocate(signsRef)
      allocate(signsRef(lenof_sign,size))
      
    end subroutine reallocate_ilutRefAdi

!------------------------------------------------------------------------------------------!

    pure function ilut_not_in_list(ilut, buffer, buffer_size) result(t_is_new)
      use DetBitOps, only: DetBitEQ
      use bit_rep_data, only: NIfD
      implicit none
      integer, intent(in) :: buffer_size
      integer(n_int), intent(in) :: ilut(0:NIfTot), buffer(0:NIfTot,buffer_size)
      logical :: t_is_new
      integer :: i_ilut

      t_is_new = .true.     
      do i_ilut = 1, buffer_size
         if(DetBitEQ(ilut, buffer(:,i_ilut), NIfD)) t_is_new = .false.
      enddo
    end function ilut_not_in_list

!------------------------------------------------------------------------------------------!

    function cpIndex(i,nRefs,cp) result(index)
      implicit none
      integer, intent(in) :: i, nRefs, cp
      integer :: index
      
      index = mod(i/(nRefs**(cp-1)),nRefs)+1
    end function cpIndex

!------------------------------------------------------------------------------------------!

    function maxProdEx(nRefs, prodLvl) result(nProdsMax)
      implicit none
      integer, intent(in) :: nRefs, prodLvl
      integer :: nProdsMax
      integer :: cLvl

      ! Gets the maximum number of possible product excitations of nRefs states up 
      ! to a product level of prodLvl
      nProdsMax = 0
      do cLvl = 2, prodLvl
         nProdsMax = nProdsMax + nRefs**cLvl
      enddo

    end function maxProdEx

!------------------------------------------------------------------------------------------!

    subroutine enable_adi
      use CalcData, only: tAllDoubsInitiators
      implicit none

      write(iout,'()') 
      write(iout,*) "Setting all double excitations to initiators"
      write(iout,*) "Using static references"
      write(iout,*) "Further notification on additional references will be given"
      write(iout,'()')
      tAllDoubsInitiators = .true.
      
      ! tSinglePartPhase = .false.
      
    end subroutine enable_adi

!------------------------------------------------------------------------------------------!

    function test_ref_double(ilut, run) result(is_accessible)
      ! We check if a target determinant for a spawn is valid in the sense
      ! that we allow any non-initiator to spawn there.
      ! This is an experimental and potentially dangerous feature as it can
      ! lead to sign instabilities
      use CalcData, only: tAccessibleDoubles, tAccessibleSingles
      use DetBitOps, only: FindBitExcitLevel
      implicit none
      integer(n_int), intent(in) :: ilut(0:NIfTot)
      integer, intent(in) :: run
      logical :: is_accessible
      integer :: iRef, exLevel
      
      is_accessible = .false.
      exLevel = -1
      do iRef = 1, nRefs
         exLevel = FindBitExcitLevel(ilutRefAdi(:,run,iRef), ilut)
         if(exLevel == 0) is_accessible = .true.
         if(tAccessibleDoubles .and. exLevel == 2) is_accessible = .true.
         if(tAccessibleSingles .and. exLevel == 1) is_accessible = .true.
         if(is_accessible) exit
      end do
    end function test_ref_double

!------------------------------------------------------------------------------------------!

    function check_sign_coherence(ilut, ilut_sign, iRef) result(is_coherent)
      use bit_reps, only: decode_bit_det
      use bit_rep_data, only: extract_sign
      use Determinants, only: get_helement
      use SystemData, only: tHPHF
      use hphf_integrals, only: hphf_off_diag_helement
      implicit none
      integer(n_int), intent(in) :: ilut(0:NIfTot)
      integer, intent(in) :: iRef
      real(dp), intent(in) :: ilut_sign(lenof_sign)
      logical :: is_coherent
      integer :: run, nI(nel), nJRef(nel)
      real(dp) :: signRef(lenof_sign)
#ifdef __CMPLX 
      complex(dp) :: tmp
#endif
      HElement_t(dp) :: h_el
      
      is_coherent = .true.
      ! Obviously, we need the sign to compare coherence
      call extract_sign(ilutRefAdi(:,1,iRef), signRef)
      
      ! For getting the matrix element, we also need the determinants
      call decode_bit_det(nI, ilut)
      call decode_bit_det(nJRef, ilutRefAdi(:,1,iRef))

      ! Then, get the matrix element
      if(tHPHF) then
         h_el = hphf_off_diag_helement(nI,nJRef(:),ilut,ilutRefAdi(:,1,iRef))
      else
         h_el = get_helement(nI,nJRef(:),ilut,ilutRefAdi(:,1,iRef))
      endif
      ! if the determinants are not coupled, ignore them
      if(abs(h_el) < eps) return

      do run = 1, inum_runs
         ! If the new ilut has sign 0, there is no need to do any check on this run
         if(mag_of_run(ilut_sign,run) > eps) then
#ifdef __CMPLX
! The complex coherence check is more effortive, so only do it in complex builds
            ! Get the relative phase of the signs
            tmp = cmplx(signRef(min_part_type(run)),signRef(max_part_type(run)),dp) / &
                 cmplx(ilut_sign(min_part_type(run)), ilut_sign(max_part_type(run)), dp)
            ! and compare it to that of H
            if(aimag(h_el*tmp) > eps .or. real(h_el*tmp) > 0.0_dp) then
               is_coherent = .false.
               return
            endif
#else
            ! For the real comparison, we just compare the signs
            if(signRef(run)*ilut_sign(run) * h_el > 0.0_dp) then
               is_coherent = .false.
               return
            endif
#endif
         endif
      enddo
    end function check_sign_coherence

!------------------------------------------------------------------------------------------!

    subroutine reset_coherence_counter()
      use FciMCData, only: nCoherentDoubles, nIncoherentDets, nCoherentSingles
      implicit none
      
      nCoherentSingles = 0
      nCoherentDoubles = 0
      nIncoherentDets = 0
    end subroutine reset_coherence_counter

!------------------------------------------------------------------------------------------!
    
    subroutine update_first_reference(ilut)
      implicit none
      integer(n_int) :: ilut(0:NIfTot)
      integer :: run
      
      do run = 1, inum_runs
         ilutRefAdi(:,run,1) = ilutRef(:,1)
      end do
      call update_single_ref_sign(1)
    end subroutine update_first_reference

!------------------------------------------------------------------------------------------!

    function giovannis_check(ilut) result(init)
      use CalcData, only: g_markers, g_markers_num
      use DetBitOps, only: CountBits
      use bit_rep_data, only: NIfD
      ! Very much like DetBitOps::FindBitExcitLevel, but it does the check for the active-space
      ! initiator criterium
      implicit none
      integer(n_int), intent(in) :: ilut(0:NIfD)
      integer(n_int) :: ilut_tmp(0:NIfD)
      integer :: nOrbs
      logical :: init

      init = .false.
      ! Get the bitwise equivalence of the input with the reference
      ilut_tmp = NOT(IEOR(ilut,ilutRef(:,1)))
      ! Now compare ilut_tmp with the markers (these are 0 for the core orbitals and 1 else)
      ilut_tmp = IAND(ilut_tmp, g_markers)
      nOrbs = CountBits(ilut_tmp,NIfD)
      ! If nOrbs is the markersize, it is an initiator
      init = (nOrbs == g_markers_num)

    end function giovannis_check
    
end module adi_references
