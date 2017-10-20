#include "macros.h"
module adi_references
use Parallel_neci
use FciMCData, only: ilutRef
use adi_data, only: ilutRefAdi, nRefs, nIRef, signsRef, nRefsSings, nRefsDoubs, &
     nTZero, SIHash, tAdiActive, tSetupSIs, NoTypeN, superInitiatorLevel, tSetupSIs, &
     tReferenceChanged, SIThreshold, tUseCaches, nIRef, signsRef, exLvlRef
use CalcData, only: InitiatorWalkNo
use bit_rep_data, only: niftot, nifdbo, extract_sign
use bit_reps, only: decode_bit_det
use DetBitOps, only: FindBitExcitLevel
use constants
use SystemData, only: nel

implicit none

interface get_sign_op
   module procedure get_sign_op_min
   module procedure get_sign_op_run
end interface get_sign_op

contains

  subroutine setup_reference_space(tPopPresent)
    use adi_data, only: tReadRefs, tProductReferences, nExProd
    use LoggingData, only: ref_filename
    implicit none
    logical, intent(in) :: tPopPresent
    integer :: nRead, nRCOld
    logical :: tGen

    if(tAdiActive) then
       ! A first guess at the number of references
       nRefs = max(nRefsDoubs, nRefsSings)

       ! If we actually did something
       tGen = .false.
       ! First, generate the reference space of nRefs determinants from the population (if present)
       if(tPopPresent) then
          call generate_ref_space()
          tGen = .true.
          if(.not. tReadRefs)&
               call print_reference_notification(1,nRefs, "Superinitiators from population")
       endif

       if(tReadRefs) then
          ! Then, add the references from the file
          call read_in_refs(ref_filename, nRead, tPopPresent)
          ! If we added references, note this (this means that the nRefs was not
          ! specified and is taken from the file)
          if(nRead > nRefs) nRefs = nRead
          tGen = .true.

          ! Print the read-in references
          call print_reference_notification(1,nRead,"Read in superinitiators")
          call print_reference_notification(nRead+1,nRefs,"Superinitiators from population")
       endif
       ! Then, add the product references
       if(tGen) then
          if(tProductReferences) then
             nRCOld = nRefs
             call add_product_references(nRefs, nExProd)
             call update_ref_signs()
             ! And prompt the output message
             call print_reference_notification(nRCOld + 1, nRefs, &
                  "Superinitiators created from excitation products")
          endif
       else
          ! If we did not do anything, only take one reference
          call reallocate_ilutRefAdi(1)
          ilutRefAdi(:,1) = ilutRef(:,1)
          nRefs = 1
          call print_reference_notification(1,1, &
               "Using only the reference determinant as superinitiator")
       endif

       call assign_SIHash_TZero()

       ! These are now the t-0 references
       nTZero = nRefs

       ! If we also used t-n (n>0) references (now called superinitiators), add them now
       if(superInitiatorLevel > 0) call add_derived_refs()

       call fill_adi_caches()

       tSetupSIs = .true.
       tReferenceChanged = .true.
    endif
   
  end subroutine setup_reference_space

!------------------------------------------------------------------------------------------!

  subroutine read_in_refs(filename, nRead, tPopPresent)
    use util_mod, only: get_free_unit
    use adi_data, only: tDelayGetRefs
    implicit none
    integer, intent(out) :: nRead
    character(255), intent(in) :: filename
    logical, intent(in) :: tPopPresent
    integer :: iunit, i, stat
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
       read(iunit, *, iostat=stat) ilutRefAdi(:,i)
       ! If there are no more dets to be read, exit
       if(stat < 0) exit
       ! If we successfully read, log it
       nRead = nRead + 1
    enddo

    close(iunit)
    ! If we read in less than nRefs detererminants, we want to get more later
    if(nRead < max(nRefsSings, nRefsDoubs) .and. .not. tPopPresent) then
       tDelayGetRefs = .true.
       nRefs = nRead
    endif
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
      implicit none
      integer(n_int), intent(inout) :: mpi_buf(0:NIfTot,nRefs)
      integer :: i, k
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
      ilutRefAdi(:,1:nRefs) = mpi_buf(:,1:nRefs)

      ! if needed, update the signs of ilutRefAdi
      if(get_sign) call update_ref_signs()
    end subroutine set_additional_references

!------------------------------------------------------------------------------------------!

    subroutine output_reference_space(filename)
      use ParallelHelper, only: root
      use util_mod, only: get_free_unit
      implicit none
      character(255), intent(in) :: filename
      integer :: iunit, i, ierr
      
      if(iProcIndex == root) then
         iunit = get_free_unit()
         open(iunit, file = filename, status = 'replace')
         ! write the references into the file
         do i = 1, nRefs
            write(iunit, *) ilutRefAdi(:,i)
         enddo
         call neci_flush(iunit)
         close(iunit)
      endif
      call MPIBarrier(ierr)
    end subroutine output_reference_space

!------------------------------------------------------------------------------------------!

    subroutine print_reference_notification(iStart, iEnd, title, legend)      
      implicit none
      integer, intent(in) :: iStart, iEnd
      character(*), intent(in) :: title
      logical, optional :: legend
      integer :: i
      
      if(iProcIndex==root .or. .true.) then
         write(iout,*) title
         if(present(legend)) write(iout,"(4A25)") &
              ! TODO: Adapt legend for multiple runs
              "Determinant (bitwise)", "Excitation level", "Coherence parameter", "Number of walkers"
         do i=iStart, iEnd
            call print_bit_rep(ilutRefAdi(:,i))
         enddo
      endif
    end subroutine print_reference_notification

!------------------------------------------------------------------------------------------!
    
    subroutine print_bit_rep(ilut)
      implicit none
      integer(n_int), intent(in) :: ilut(0:NIfTot)
      real(dp) :: temp_sgn(lenof_sign)
      integer :: j
      ! First, print the determinant (bitwise)
      call WriteDetBit(iout,ilut,.false.)
      ! Then, the excitation level
      write(iout, "(5X)", advance = 'no')
      write(iout, "(G1.4)", advance = 'no') FindBitExcitLevel(ilut, ilutRef(:,1))
      ! And the sign coherence parameter
      write(iout, "(G16.7)", advance = 'no') get_sign_op_min(ilut)
      ! And then the sign
      call extract_sign(ilut,temp_sgn)
      do j = 1, lenof_sign
         write(iout,"(G16.7)",advance='no') temp_sgn(j)
      enddo
      
      write(iout,'()')

    end subroutine print_bit_rep

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
      use hash, only: hash_table_lookup
      use FciMCData, only: CurrentDets, HashIndex
      use Parallel_neci, only: MPISumAll
      implicit none
      integer, intent(in) :: iRef
      integer:: nI(nel), hash_val, index
      real(dp) :: tmp_sgn(lenof_sign), mpi_sgn(lenof_sign)
      logical :: tSuccess      
      
      call decode_bit_det(nI,ilutRefAdi(:,iRef))
      call hash_table_lookup(nI, ilutRefAdi(:,iRef),NIfDBO,HashIndex,CurrentDets,&
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
      call encode_sign(ilutRefAdi(:,iRef),mpi_sgn)

    end subroutine update_single_ref_sign

!------------------------------------------------------------------------------------------!

    subroutine spin_symmetrize_references()
      use DetBitOps, only: spin_flip, DetBitEQ
      ! If using HPHF, we want to have the spinflipped version of each reference, too
      ! This is a bit cumbersome to implement, as the spinflipped version might or
      ! might not be already in the references
      implicit none
      integer :: iRef, jRef,  cRef
      integer(n_int) :: ilutRef_new(0:NIfTot,2*nRefs), tmp_ilut(0:NIfTot)
      logical :: missing
      
      cRef = 1
      do iRef = 1 ,nRefs
         ! First, copy the current ilut to a new buffer
         ilutRef_new(:,cRef) = ilutRefAdi(:,iRef)
         cRef = cRef + 1
         ! get the spinflipped determinant
         tmp_ilut = spin_flip(ilutRefAdi(:,iRef))
         missing = .true.
         ! check if it is already present
         do jRef = 1, nRefs
            if(DetBitEQ(tmp_ilut,ilutRefAdi(:,jRef),NIfDBO)) missing = .false.
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
      ilutRefAdi(:,1:cRef) = ilutRef_new(:,1:cRef)
    end subroutine spin_symmetrize_references

!------------------------------------------------------------------------------------------!

    subroutine add_product_references(nRefsIn, prodLvl)
      implicit none
      integer, intent(in) :: nRefsIn
      integer, intent(in) :: prodLvl
      integer(n_int) :: tmp_ilut(0:NIfTot)
      integer(n_int), allocatable :: prod_buffer(:,:)
      integer :: nProdsMax, nProds, i, cLvl
      logical :: t_is_valid

      ! Get the maximum number of product excitations
      nProdsMax = maxProdEx(nRefsIn, prodLvl)

      ! Setup the buffer
      allocate(prod_buffer(0:NifTot,nProdsMax))
      prod_buffer(:,1:nRefsIn) = ilutRefAdi(:,1:nRefsIn)
      nProds = nRefsIn

      do cLvl = 2, prodLvl
         ! For each product level, get all possible product excitations
         do i = 1, nRefsIn**cLvl
            ! Get the excitation product corresponding to i
            call get_product_excitation(i, cLvl, nRefsIn, tmp_ilut, t_is_valid)
            
            ! Check, if it is already in the list (e.g. if both tau operators are those of
            ! the reference)
            if(t_is_valid) t_is_valid = ilut_not_in_list(tmp_ilut, prod_buffer, nProds)
         
            ! If all excitations were valid and the new ilut is not already in the buffer, 
            ! add it
            if(t_is_valid) then
               nProds = nProds + 1
               prod_buffer(:,nProds) = tmp_ilut
            endif
         end do
         ! Do this for each targeted product level
      end do
            
      ! Write the buffer to ilutRefAdi
      call reallocate_ilutRefAdi(nProds)
      ! Reassign the number of references
      nRefs = nProds
      ilutRefAdi = prod_buffer(:,1:nProds)

      deallocate(prod_buffer)
    end subroutine add_product_references

!------------------------------------------------------------------------------------------!

    subroutine get_product_excitation(i, cLvl, nRefs, tmp_ilut, isValid)
      use DetBitOps, only: CountBits
      use bit_rep_data, only: NIfD
      implicit none
      integer, intent(in) :: i, cLvl, nRefs
      integer(n_int), intent(out) :: tmp_ilut(0:NIfTot)
      logical, intent(out) :: isValid
      integer(n_int) :: tau(0:NIfTot), tau_cc(0:NIfTot)
      integer :: cp

      tmp_ilut = ilutRef(:,1)
      isValid = .true.
      tau_cc = 0_n_int
      ! By getting the cLvl excitation operators that correspond to i
      do cp = 1, cLvl
         ! Store the excitation operators in tau
         tau = IEOR(ilutRefAdi(:,cpIndex(i,nRefs,cp)),ilutRef(:,1))

         ! Check if we do not annihilate something twice
         if(any(IAND(tau_cc(0:NIfD),tau(0:NIfD)) .ne. 0_n_int)) then
            isValid = .false.
            exit
         endif
         ! And apply it on the temporary iLut
         tau_cc(0:NIfD) = IEOR(tau_cc(0:NIfD),tau(0:NIfD))

      end do

      ! Check if it's valid, if not, go to the next value of i
      if(CountBits(tmp_ilut,NIfD) .ne. nEl) isValid = .false.
      tmp_ilut(0:NIfD) = IEOR(tmp_ilut(0:NIfD), tau_cc(0:NIfD))
    end subroutine get_product_excitation

!------------------------------------------------------------------------------------------!

    subroutine add_derived_refs()
      implicit none
      integer :: level, nRCOld
      character :: type

      ! Add sign-coherent states as superinitiators
      do level = 1, superInitiatorLevel
         nRCOld = nRefs
         call generate_type_n_refs()
         write(type,'(i1)') level
         ! Print the type-n superinitiators of the current level
         call print_reference_notification(nRCOld + 1, nRefs, "Type-"//type//&
              " superinitiators")
      enddo
    end subroutine add_derived_refs

!------------------------------------------------------------------------------------------!

    subroutine generate_type_n_refs()
      use FciMCData, only: CurrentDets, TotWalkers
      implicit none
      integer :: i, nTOne, nBlocks, ierr, run
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
         ! Check if a determinant should be a type-1 reference (on any run)
         is_tone = .false.
         do run = 1, inum_runs
            is_tone = is_tone .or. check_type_n_ref(CurrentDets(:,i), tmp_sign, run)
            ! TODO: Set the flag for the corresponding run
         enddo
         if(is_tone) then
            nTOne = nTOne + 1
            ! If we exceed the memory of the buffer, add another block
            if(nTOne > nBlocks*allocBlock) then
               nBlocks = nBlocks + 1
               call resize_ilut_list(refBuf, (nBlocks-1) * allocBlock ,nBlocks * allocBlock)
               if(ierr .ne. 0) call stop_all(this_routine, "Unable to allocate temporary list")
            endif
            ! Add the determinant to the temporary list
            refBuf(:,nTOne) = CurrentDets(:,i)
         endif
      end do

      ! And construct the array of new superinitiators
      call add_type_n_refs(refBuf, nTOne)

      deallocate(refBuf)
      
    end subroutine generate_type_n_refs

!------------------------------------------------------------------------------------------!

    subroutine add_type_n_refs(list, listSize)
      implicit none
      integer, intent(in) :: listSize
      integer(n_int), intent(in) :: list(0:NIfTot,listSize)
      integer(n_int), allocatable :: mpi_buf(:,:), buffer(:,:)
      integer :: nRCOld, i, nNew, all_refs_found, ierr
      integer(MPIArg) :: refs_found_per_proc(0:nProcessors-1), refs_displs(0:nProcessors-1)
      logical :: tSuccess
      integer(MPIArg) :: mpi_refs_found
      
      ! First, communicate the list between the processors
      mpi_refs_found = int(listSize,MPIArg)
      call MPIAllGather(mpi_refs_found, refs_found_per_proc, ierr)
      all_refs_found = sum(refs_found_per_proc)

      ! We only now know the required size of the temporaries
      allocate(mpi_buf(0:NIfTot, all_refs_found))
      allocate(buffer(0:NIfTot, all_refs_found))

      refs_displs(0) = 0
      do i = 1, nProcessors - 1
         refs_displs(i) = sum(refs_found_per_proc(0:i-1))
      enddo
      ! Store them on all processors
      call MPIAllGatherV(list(0:NIfTot, 1:listSize), mpi_buf, &
           refs_found_per_proc, refs_displs)

      nRCOld = nRefs
      nNew = 0
      ! we need to allocate the buffer to maximum size because adding to the
      ! hashtable requires the target array (here: ilutRefAdi) to have sufficient size
      call resize_ilutRefAdi(nRefs + all_refs_found)
      ! note that we only need the allocation, the reassignment of nRefs is pointless
      ! First, pick those iluts from list, which are not already present in ilutRefAdi
      do i = 1, all_refs_found
         call add_superinitiator_to_hashtable(mpi_buf(:,i), nRCOld + nNew + 1, tSuccess)
         if(.not. tSuccess) then
            nNew = nNew + 1
            buffer(:,nNew) = mpi_buf(:,i)
         endif
      enddo
      ! resize the ilutRefAdi array
      call resize_ilutRefAdi(nRCOld + nNew)
      ! and add the new iluts
      ilutRefAdi(:,(nRCOld + 1):(nRCOld + nNew)) = buffer(:,1:nNew)

      ! and then check the self-consistency condition, keeping the old SIs
      if(nNew  >0) call type_n_self_consistency_loop(nRCOld)

      deallocate(mpi_buf)
      deallocate(buffer)
    end subroutine add_type_n_refs

!------------------------------------------------------------------------------------------!

    function check_type_n_ref(ilut, ilut_sign, run) result(is_tone)
      implicit none
      integer(n_int), intent(in) :: ilut(0:NIfTot)
      ! We pass the sign as an extra argument to prevent redundant extraction
      real(dp), intent(in) :: ilut_sign(lenof_sign)
      integer, intent(in) :: run
      real(dp) :: xi
      logical :: is_coherent
      logical :: is_tone
      
      is_tone = .false.
      is_coherent = .false.

      ! Only consider those determinants with a sufficient population
      if(mag_of_run(ilut_sign, run) < NoTypeN) return

      ! obtain the xi-parameter
      xi = get_sign_op_run(ilut, run) 
      ! and compare it to the superinitiator threshold
      if(xi > SIThreshold) is_tone = .true.
    end function check_type_n_ref
    
!------------------------------------------------------------------------------------------!

    subroutine type_n_self_consistency_loop(nKeep)
      implicit none
      integer, intent(in) :: nKeep
      integer :: i

      ! we can remove at most nrefs SIs
      do i = 1, nRefs
         ! once we do not remove any SIs, self-consistency is reached
         if(type_n_self_consistency_step(nKeep)) exit         
      enddo

      ! We need to setup the hashtable anew because the SIs are reordered in
      ! the process of self-consistency loop
      call assign_SIHash_TZero()
    end subroutine type_n_self_consistency_loop

!------------------------------------------------------------------------------------------!
    
    function type_n_self_consistency_step(nKeep) result(valid)
      implicit none
      integer, intent(in) :: nKeep
      logical :: valid
      integer :: iRef, min_iRef
      real(dp) :: sub_xi, min_xi
      
      valid = .true.
      min_xi = 1.0_dp
      min_iRef = nKeep + 1

      ! The first nKeep references are not touchable         
      do iRef = nKeep + 1, nRefs
         ! Get the minimal xi
         sub_xi = get_sign_op_min(ilutRefAdi(:,iRef))
         ! And the corresponding SI
         if(sub_xi < min_xi) then
            min_xi = sub_xi
            min_iRef = iRef
         endif
      end do

      ! But we need to throw something out if the first nKeep SIs fall below threshold
      do iRef = 1, nKeep
         if(get_sign_op_min(ilutRefAdi(:,iRef)) < SIThreshold) valid = .false.
      enddo
      
      ! and throw out the SI with minimal xi if below threshold or if there is a 
      ! threshold violation in the first nKeep entries
      if(min_xi < SIThreshold .or. (.not. valid)) then
         call remove_superinitiator(min_iRef)
         valid = .false.
      endif
      
    end function type_n_self_consistency_step

!------------------------------------------------------------------------------------------!

    function check_sign_coherence(ilut, nI, ilut_sign, iRef, run) result(is_coherent)
      use Determinants, only: get_helement
      use SystemData, only: tHPHF
      use hphf_integrals, only: hphf_off_diag_helement
      implicit none
      integer(n_int), intent(in) :: ilut(0:NIfTot)
      integer, intent(in) :: iRef, run, nI(nel)
      real(dp), intent(in) :: ilut_sign(lenof_sign)
      logical :: is_coherent
      integer :: nJRef(nel)
      real(dp) :: signRef(lenof_sign)
#ifdef __CMPLX 
      complex(dp) :: tmp
#endif
      HElement_t(dp) :: h_el
      
      is_coherent = .true.
      ! Obviously, we need the sign to compare coherence
      call extract_sign(ilutRefAdi(:,iRef), signRef)
      
      ! For getting the matrix element, we also need the determinants
      call decode_bit_det(nJRef, ilutRefAdi(:,iRef))

      ! Then, get the matrix element
      if(tHPHF) then
         h_el = hphf_off_diag_helement(nI,nJRef(:),ilut,ilutRefAdi(:,iRef))
      else
         h_el = get_helement(nI,nJRef(:),ilut,ilutRefAdi(:,iRef))
      endif
      ! if the determinants are not coupled, ignore them
      if(abs(h_el) < eps) return

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
    end function check_sign_coherence

!------------------------------------------------------------------------------------------!
  
  subroutine initialize_c_caches(signedCache, unsignedCache)
    implicit none
    HElement_t(dp), intent(out) :: signedCache
    real(dp), intent(out) :: unsignedCache
    
    ! We potentially want to add the diagonal terms here
    signedCache = 0.0_dp
    unsignedCache = 0.0_dp
  end subroutine initialize_c_caches

!------------------------------------------------------------------------------------------!

  subroutine update_coherence_check(ilut, nI, i, run, signedCache, unsignedCache)
    use SystemData, only: tHPHF 
    use Determinants, only: get_helement
    use hphf_integrals, only: hphf_off_diag_helement
    implicit none
    integer, intent(in) :: nI(nel), i, run
    integer(n_int), intent(in) :: ilut(0:NIfTot)
    HElement_t(dp), intent(inout) :: signedCache
    real(dp), intent(inout) :: unsignedCache
    HElement_t(dp) :: h_el, tmp
    character(*), parameter :: this_routine = "upadte_coherence_check"

    ! TODO: Only if ilutRefAdi(:,i) is a SI on this run

    ! We require cached data here
    if(.not. tUseCaches) call fill_adi_caches
    ! First, get the matrix element
    if(tHPHF) then
       h_el = hphf_off_diag_helement(nI,nIRef(:,i),ilut,ilutRefAdi(:,i))
    else
       h_el = get_helement(nI,nIRef(:,i),ilut,ilutRefAdi(:,i))
    endif
    ! Only proceed if the determinants are coupled
    if(abs(h_el) < eps) return

    ! Add tmp = Hij cj to the caches
#ifdef __CMPLX
    tmp = h_el * cmplx(signsRef(min_part_type(run),i),&
         signsRef(max_part_type(run),i),dp)
#else
    tmp = h_el * signsRef(run,i)
#endif
    signedCache = signedCache + tmp
    unsignedCache = unsignedCache + abs(tmp)
  end subroutine update_coherence_check

!------------------------------------------------------------------------------------------!

  subroutine eval_coherence(signedCache, unsignedCache, sgn, staticInit)
    ! Note that tweakcoherentdoubles and tavcoherentdoubles are mutually exclusive 
    ! in the current implementation
    use adi_data, only: tWeakCoherentDoubles, tAvCoherentDoubles, coherenceThreshold
    implicit none
    HElement_t(dp), intent(in) :: signedCache
    real(dp), intent(in) :: unsignedCache, sgn
    logical, intent(inout) :: staticInit

    ! Only need to check if we are looking at a double
    if(unsignedCache > eps) then
       ! We disable superinitiator-related initiators if they fail the coherence check
       ! else, we leave it as it is
       if(tWeakCoherentDoubles) then
          if(abs(signedCache) < coherenceThreshold*unsignedCache) staticInit = .false.
       endif

       ! If we do averaged coherence check, we check versus the sign of the ilut
       ! We recommend using both, av and weak
       if(tAvCoherentDoubles) then
          if(real(signedCache * sgn,dp) > 0.0_dp) staticInit = .false.
       endif
    endif

  end subroutine eval_coherence

!------------------------------------------------------------------------------------------!

  function get_sign_op_min(ilut) result(xi)
    implicit none
    integer(n_int) :: ilut(0:NIfTot)
    integer :: run
    real(dp) :: min_xi, xi

    min_xi = 1.0_dp
    do run = 1, inum_runs
       xi = get_sign_op(ilut, run)
       if(xi < min_xi) min_xi = xi
    enddo
  end function get_sign_op_min

!------------------------------------------------------------------------------------------!

  function get_sign_op_run(ilut, run) result(xi)
    ! compute the sign problem indicator xi for a given determinant
    implicit none
    integer(n_int), intent(in) :: ilut(0:NIfTot)
    integer, intent(in) :: run
    real(dp) :: xi
    integer :: iRef, nI(nel), exLevel
    real(dp) :: unsignedCache
    HElement_t(dp) :: signedCache
    
    call initialize_c_caches(signedCache, unsignedCache)
    ! Sum up all Hij cj for all superinitiators j
    do iRef = 1, nRefs
       exLevel = FindBitExcitLevel(ilutRefAdi(:,iRef), ilut)
       ! Of course, only singles/doubles of ilut can contribute
       if(exLevel < 3) then
          call decode_bit_det(nI, ilut)
          call update_coherence_check(ilut, nI, iRef, run, &
               signedCache, unsignedCache)
       endif
    enddo
    ! Get the ratio of the caches, which is xi
    if(abs(unsignedCache) > eps) then
       xi = abs(signedCache) / unsignedCache
    else
       xi = 0.0_dp
    endif
  end function get_sign_op_run

!------------------------------------------------------------------------------------------!

  subroutine fill_adi_caches()
    implicit none
    integer :: iRef

    call allocate_adi_caches()
    
    do iRef = 1, nRefs
       call decode_bit_det(nIRef(:,iRef),ilutRefAdi(:,iRef))
       call extract_sign(ilutRefAdi(:,iRef), signsRef(:,iRef))
       exLvlRef(iRef) = FindBitExcitLevel(ilutRef(:,1), ilutRefAdi(:,iRef))
    end do

    tUseCaches = .true.
  end subroutine fill_adi_caches
!------------------------------------------------------------------------------------------!

  subroutine allocate_adi_caches()
    implicit none

    call deallocate_adi_caches()

    allocate(nIRef(nel, nRefs))
    allocate(signsRef(lenof_sign, nRefs))
    allocate(exLvlRef(nRefs))
  end subroutine allocate_adi_caches

!------------------------------------------------------------------------------------------!

  subroutine deallocate_adi_caches()
    implicit none
    
    if(allocated(nIRef)) deallocate(nIRef)
    if(allocated(signsRef)) deallocate(signsRef)
    if(allocated(exLvlRef)) deallocate(exLvlRef)
  end subroutine deallocate_adi_caches

!------------------------------------------------------------------------------------------!

    subroutine resize_ilutRefAdi(newSize)
      implicit none
      integer, intent(in) :: newSize

      ! see resize_ilut_list for the implementation
      call resize_ilut_list(ilutRefAdi, nRefs, newSize)      
      ! The new number of references
      nRefs = newSize

      tUseCaches = .false.
         
    end subroutine resize_ilutRefAdi

!------------------------------------------------------------------------------------------!

    subroutine resize_ilut_list(list, oldSize, newSize)
      implicit none
      integer, intent(in) :: newSize, oldSize
      integer(n_int), allocatable, intent(inout) :: list(:,:)
      integer(n_int) :: tmp(0:NIfTot,oldSize)
      integer :: ierr
      logical :: tUseTmp
      character(*), parameter :: this_routine = "resize_ilut_list"

      ! We store the current list in a temporary, if existent
      tUseTmp = .false.
      ! If the list is already allocated, clear it
      if(allocated(list)) then
         tmp(:,:) = list(:,:)
         tUseTmp = .true.
         deallocate(list)
      endif

      allocate(list(0:NIfTot,newSize), stat = ierr)
      if(ierr .ne. 0) call stop_all(this_routine, "Not enough memory to resize ilut list")
      list = 0
      
      ! Write the old entries of list into the new memory
      if(oldSize < newSize) then
         list(:,1:oldSize) = tmp(:,1:oldSize)
      else
         list(:,1:newSize) = tmp(:,1:newSize)
      endif
      
    end subroutine resize_ilut_list

!------------------------------------------------------------------------------------------!

    subroutine reallocate_ilutRefAdi(size)
      implicit none
      integer, intent(in) :: size
      
      ! reallocate first the ilutref itself
      if(allocated(ilutRefAdi)) deallocate(ilutRefAdi)
      allocate(ilutRefAdi(0:NIfTot,size))
      ilutRefAdi = 0_n_int
      
      tSetupSIs = .false.
      tUseCaches = .false.
      
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
      use adi_data, only: tAllDoubsInitiators, tAdiActive
      implicit none

      write(iout,'()') 
      write(iout,*) "Setting all double excitations to initiators"
      write(iout,*) "Using static references"
      write(iout,*) "Further notification on additional references will be given"
      write(iout,'()')
      tAllDoubsInitiators = .true.
      tAdiActive = .true.
      
      ! tSinglePartPhase = .false.
      
    end subroutine enable_adi

!------------------------------------------------------------------------------------------!

    function test_ref_double(ilut, run) result(is_accessible)
      ! We check if a target determinant for a spawn is valid in the sense
      ! that we allow any non-initiator to spawn there.
      ! This is an experimental and potentially dangerous feature as it can
      ! lead to sign instabilities
      use adi_data, only: tAccessibleDoubles, tAccessibleSingles
      implicit none
      integer(n_int), intent(in) :: ilut(0:NIfTot)
      integer, intent(in) :: run
      logical :: is_accessible
      integer :: iRef, exLevel
      
      is_accessible = .false.
      exLevel = -1
      do iRef = 1, nRefs
         exLevel = FindBitExcitLevel(ilutRefAdi(:,iRef), ilut)
         if(exLevel == 0) is_accessible = .true.
         if(tAccessibleDoubles .and. exLevel == 2) is_accessible = .true.
         if(tAccessibleSingles .and. exLevel == 1) is_accessible = .true.
         if(is_accessible) exit
      end do
    end function test_ref_double


!------------------------------------------------------------------------------------------!

    subroutine reset_coherence_counter()
      use adi_data, only: nCoherentDoubles, nIncoherentDets, nCoherentSingles
      implicit none
      
      nCoherentSingles = 0
      nCoherentDoubles = 0
      nIncoherentDets = 0
    end subroutine reset_coherence_counter

!------------------------------------------------------------------------------------------!
    
    subroutine update_first_reference()
      use FciMCData, only: tSinglePartPhase
      implicit none

      call setup_reference_space(all(.not. tSinglePartPhase))
    end subroutine update_first_reference

!------------------------------------------------------------------------------------------!

    subroutine setup_SIHash()
      use adi_data, only:  htBlock
      use FciMCData, only: HashLengthFrac
      use hash, only: init_hash_table
      implicit none
      
      htBlock = 5*HashLengthFrac*nRefs
      if(associated(SIHash)) deallocate(SIHash)
      allocate(SIHash(htBlock))
      call init_hash_table(SIHash)
    end subroutine setup_SIHash

!------------------------------------------------------------------------------------------!

    subroutine remove_superinitiator(iRef)
      ! calling this once changes the indices of the entries, the hash-table
      ! is hence invalidated
      implicit none
      integer, intent(in) :: iRef
      integer(n_int) :: tmp(0:NIfTot)

      ! then, shrink ilutRefAdi
      ! as the last entry will be deleted, move iRef to the back
      if(iRef <  nRefs) then
         tmp = ilutRefAdi(:,nRefs)
         ilutRefAdi(:,nRefs) = ilutRefAdi(:,iRef)
         ilutRefAdi(:,iRef) = tmp
      endif
      call resize_ilutRefAdi(nRefs - 1)

    end subroutine remove_superinitiator

!------------------------------------------------------------------------------------------!

    subroutine add_superinitiator_to_hashtable(ilut, cRef, tSuccess)
      use hash, only: add_hash_table_entry, hash_table_lookup
      use bit_reps, only: decode_bit_det
      implicit none
      integer(n_int), intent(in) :: ilut(0:NIfTot)
      integer, intent(in) :: cRef
      integer :: nI(nel), index, hashVal
      logical, intent(out) :: tSuccess
      
      call decode_bit_det(nI, ilut)
      call hash_table_lookup(nI, ilut, NIfDBO, SIHash, ilutRefAdi, index, hashVal, tSuccess)
      ! If the SI is already present, do nothing
      ! TODO: Set flags for multi-replica
      if(tSuccess) return
      
      index = cRef
      ! Else, add a hashtable entry
      call add_hash_table_entry(SIHash, index, hashVal)
    end subroutine add_superinitiator_to_hashtable

!------------------------------------------------------------------------------------------!

    subroutine assign_SIHash_TZero()
      use hash, only: clear_hash_table, fill_in_hash_table
      implicit none
      integer :: iRef
      logical :: tSuccess
      character(*), parameter :: this_routine = "assign_SIHash_TZero"
      
      tSuccess = .false.
      ! make sure the hash table has an appropiate size
      call setup_SIHash()
      do iRef = 1, nRefs
         call add_superinitiator_to_hashtable(ilutRefAdi(:,iRef), iRef, tSuccess)
         ! By construction, there can't be duplicates within the type-0 SIs
         if(tSuccess) call stop_all(this_routine, "Duplicate type-0 superinitiator")
      enddo
            
    end subroutine assign_SIHash_TZero

!------------------------------------------------------------------------------------------!
    
    subroutine check_superinitiator(ilut, nI, staticInit) 
      use hash, only: hash_table_lookup
      implicit none
      integer(n_int), intent(in) :: ilut(0:NIfTot)
      integer, intent(in) :: nI(nel)
      logical, intent(inout) :: staticInit
      integer :: index, hashVal
      logical :: tSuccess
      
      ! Only perform the check if the superinitiators are initialized
      if(tSetupSIs) then
         ! Now, look into the SI-hashtable and return whether nI is there
         call hash_table_lookup(nI, ilut, NIfDBO, SIHash, ilutRefAdi, &
              index, hashVal, tSuccess)
         if(tSuccess) staticInit = .true.
      endif
    end subroutine check_superinitiator

!------------------------------------------------------------------------------------------!

    function giovannis_check(ilut) result(init)
      use adi_data, only: g_markers, g_markers_num
      use DetBitOps, only: CountBits
      use bit_rep_data, only: NIfD
      ! Very much like DetBitOps::FindBitExcitLevel, but it does the check for the active-space
      ! initiator criterium
      implicit none
      integer(n_int), intent(in) :: ilut(0:NIfD)
      integer(n_int) :: tmp_ilut(0:NIfD)
      integer :: nOrbs
      logical :: init

      init = .false.
      ! Get the bitwise equivalence of the input with the reference
      tmp_ilut = NOT(IEOR(ilut,ilutRef(:,1)))
      ! Now compare tmp_ilut with the markers (these are 0 for the core orbitals and 1 else)
      tmp_ilut = IAND(tmp_ilut, g_markers)
      nOrbs = CountBits(tmp_ilut,NIfD)
      ! If nOrbs is the markersize, it is an initiator
      init = (nOrbs == g_markers_num)

    end function giovannis_check

!------------------------------------------------------------------------------------------!

    subroutine clean_adi()
      use adi_data, only: g_markers
      implicit none
      
      call deallocate_adi_caches()
      if(associated(SIHash)) deallocate(SIHash)
      if(allocated(ilutRefAdi)) deallocate(ilutRefAdi)
      if(allocated(g_markers)) deallocate(g_markers)
    end subroutine clean_adi
    
end module adi_references
