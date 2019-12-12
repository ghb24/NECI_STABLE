#include "macros.h"
MODULE Determinants
    use constants, only: dp, n_int, bits_n_int, int64, maxExcit
    use SystemData, only: BasisFN, tCSF, nel, G1, Brr, ECore, ALat, NMSH, &
                          nBasis, nBasisMax, tStoreAsExcitations, tHPHFInts, &
                          tCSF, tCPMD, tPickVirtUniform, LMS, modk_offdiag, &
                          t_lattice_model, arr, lms, tFixLz, tUEGSpecifyMomentum, &
                          tRef_Not_HF, tMolpro, tHub, tUEG, &
                          nClosedOrbs, nOccOrbs, nIrreps, tspn, irrepOrbOffset
    use IntegralsData, only: UMat, FCK, NMAX
    use csf, only: det_to_random_csf, iscsf, csf_orbital_mask, &
                   csf_yama_bit, CSFGetHelement
    use sltcnd_mod, only: sltcnd, sltcnd_excit_old, sltcnd_compat, &
                          sltcnd_knowIC, SumFock, CalcFockOrbEnergy
    use procedure_pointers, only: sltcnd_2
    use global_utilities
    use sort_mod
    use DetBitOps, only: EncodeBitDet, count_open_orbs, spatial_bit_det
    use DeterminantData
    use bit_reps
    use MemoryManager, only: TagIntType
    use lattice_mod, only: get_helement_lattice
    use util_mod, only: NECI_ICOPY
    use SymData , only : nSymLabels,SymLabelList,SymLabelCounts,TwoCycleSymGens
    use sym_mod
    use sort_mod
    use global_utilities

    implicit none

    ! TODO: Add an interface for getting a diagonal helement with an ordered
    !       list, or with only a bit-det
    interface get_helement
        module procedure get_helement_compat
        module procedure get_helement_excit
        module procedure get_helement_normal
    end interface

    save
! Set by Calc on input
      INTEGER nActiveSpace(2)
      INTEGER, DIMENSION(:), POINTER :: SPECDET => null()
      INTEGER(TagIntType) :: tagSPECDET=0
      Logical TSPECDET

!nActiveBasis(1) is the lowest non-active orbital
!nActiveBasis(2) is the highest active orbital.  There can be virtuals above this.
!  Active orbitals are used for generating the determinants whose energy/weight is to be found
      Integer nActiveBasis(2)
!  Set by input to indicate which type of active basis we need
      Integer iActiveBasis
!Not really Locals - needed for the DetCalc
      TYPE(BasisFN) ISym
!Used to be from uhfdet.inc
      INTEGER nUHFDet(5000)
      real(dp)  E0HFDet

      INTEGER, allocatable :: DefDet(:)
      Logical :: tDefineDet
      integer(TagIntType) :: tagDefDet=0

contains

  Subroutine DetPreFreezeInit()
    Use global_utilities
    use SystemData, only : nEl, ECore, Arr, Brr, G1, nBasis, LMS, nBasisMax,&
         tFixLz, tUEGSpecifyMomentum, tRef_Not_HF
    use SystemData, only : tMolpro
    use sym_mod
    use util_mod, only: NECI_ICOPY
    use sltcnd_mod, only: CalcFockOrbEnergy
    integer ierr, ms, iEl, flagAlpha, iIrrep, msTmp
    integer i,j,Lz,OrbOrder(8,2),FDetTemp(NEl),lmsMax
    type(BasisFn) s
    logical :: tGenFDet
    HElement_t(dp) :: OrbE
    character(25), parameter :: this_routine='DetPreFreezeInit'
    Allocate(FDet(nEl), stat=ierr)
    LogAlloc(ierr, 'FDet', nEl, 4, tagFDet)
    IF(tDefineDet) THEN
       WRITE(6,*) 'Defining FDet according to input'
       do i=1,NEl
          FDet(i)=DefDet(i)
       enddo
       call assignOccOrbs()

       ! A quick check that we have defined a reasonable det.
       ms = sum(get_spin_pn(fdet(1:nel)))
       if (abs(ms) /= abs(lms) .and. .not. tCSF) then
          write(6,*) 'LMS', lms
          write(6,*) 'Calculated Ms', ms
          call stop_all (this_routine, "Defined determinant has the &
               &wrong Ms value. Change DEFINEDET or &
               &SPIN-RESTRICT")
       end if
       tRef_Not_HF = .true.
    else
       if((sum(nOccOrbs) + sum(nClosedOrbs)) .eq. nel) then
          tGenFDet = .false.
          iEl = 1
          msTmp = -1*lms
          do i = 1, nIrreps
             ! doubly occupy the closed orbs
             do j = 1, nClosedOrbs(i)
                FDet(iEl) = irrepOrbOffset(i) + 2*j - 1
                iEl = iEl + 1
                FDet(iEl) = irrepOrbOffset(i) + 2*j
                iEl = iEl + 1
             end do
             ! now distribute electrons to the open orbs
             ! such that the total sz matches the requested
             ! only consider occ orbs which are not closed
             do j = nClosedOrbs(i)+1, nOccOrbs(i)
                if(msTmp<0) then
                   flagAlpha = 0
                else
                   flagAlpha = 1
                endif
                FDet(iEl) = irrepOrbOffset(i) + 2*j - flagAlpha
                iEl = iEl + 1
                msTmp = msTmp + (1-2*flagAlpha)
             end do
          end do
          call sort(FDet)
       else
          CALL GENFDET(FDET)
          call assignOccOrbs()
          IF(tUEGSpecifyMomentum) THEN
             WRITE(6,*) 'Defining FDet according to a momentum input'
             CALL ModifyMomentum(FDET)
          ENDIF
          tRef_Not_HF = .false.
       endif
    ENDIF
    !      ENDIF
    WRITE(6,"(A)",advance='no') " Fermi det (D0):"
    call write_det (6, FDET, .true.)
    Call GetSym(FDet,nEl,G1,nBasisMax,s)
    WRITE(6,"(A)",advance='no') " Symmetry: "
    Call WriteSym(6,s%Sym,.true.)
    IF(tFixLz) THEN
       Call GetLz(FDet,nEl,Lz)
       WRITE(6,"(A,I5)") "Lz of Fermi det:",Lz
    ENDIF
    CALL NECI_ICOPY(NEL,FDET,1,NUHFDET,1)
    if(tMolpro) then
       !Orbitals are ordered by occupation number from MOLPRO, and not reordered in NECI
       !Therefore, we know HF determinant is first four occupied orbitals.
       write(6,"(A)") "NECI called from MOLPRO, so assuming orbitals ordered by occupation number."
       if(.not.tDefineDet) then
          FDetTemp(:)=FDet(:)
       else
          !We have defined our own reference determinant, but still use the first orbitals for the calculation
          !of 'orbital energies'
          CALL GENFDET(FDETTEMP)
       endif
       write(6,"(A)") "Calculating orbital energies..."
       do i=1,nBasis
          OrbE=CalcFockOrbEnergy(i,FDetTemp)
          Arr(i,1)=real(OrbE,dp)
          Brr(i)=i
       enddo
       write(6,"(A)") "Reordering basis by orbital energies..."
       OrbOrder(:,:)=0
       call ORDERBASIS(NBASIS,Arr,Brr,OrbOrder,nBasisMax,G1)
       !However, we reorder them here
       call writebasis(6,G1,nBasis,Arr,Brr)
    endif
    E0HFDET=ECORE
    DO I=1,NEL
       E0HFDET=E0HFDET+ARR(NUHFDET(i),2)
    ENDDO
    WRITE(6,*) "Fock operator energy:",E0HFDET

    ! Store the value of Ms for use in other areas
    calculated_ms = sum(get_spin_pn(fdet(1:nel)))

  contains

    subroutine assignOccOrbs
      implicit none
      integer :: k
      ! assign the occ/closed orbs
      nOccOrbs = 0
      nClosedOrbs = 0
      do k = 1, nel
         if(k>1) then
            if(is_alpha(FDet(k)) .and. FDet(k-1).eq.FDet(k)-1) then
               ! we do not need to resolve the irrep anymore (not relevant
               ! for further usage)
               nClosedOrbs(1) = nClosedOrbs(1) + 1
               cycle
            endif
         end if
         nOccOrbs(1) = nOccOrbs(1) + 1
      end do
    end subroutine assignOccOrbs
  End Subroutine DetPreFreezeInit

    Subroutine DetInit()
      real(dp) DNDET
      integer i,j
      integer(int64) nDet
      integer :: alpha,beta,symalpha,symbeta,endsymstate
      LOGICAL :: tSuccess,tFoundOrbs(nBasis)
      integer :: ncsf

      WRITE(6,*) "SYMMETRY MULTIPLICATION TABLE"
      CALL WRITESYMTABLE(6)

      CALL GENSymStatePairs(NBASIS/2,.false.)


!iActiveBasis is a copy of nPaths
      IF(iActiveBasis.eq.-2) then
!  PATHS ACTIVE SETS
         Call GenActiveBasis(ARR,nBasis,nEl,nActiveBasis,nActiveSpace(1),nActiveSpace(2))
      elseif(iActiveBasis.eq.-3) then
!  PATHS ACTIVE ORBITALS
         nActiveBasis(1)=nEl+1-nActiveSpace(1)
         nActiveBasis(2)=nEl+nActiveSpace(2)
         WRITE(6,*) "Active space:", nActiveBasis(1)," TO ", nActiveBasis(2)," (ordered labels)."
         WRITE(6,*) "Active space electrons:",nEl-nActiveBasis(1)+1
      else
         nActiveBasis(1)=1
         nActiveBasis(2)=nBasis
      endif
!C.. Work out a preliminary Fermi det
!      IF(FDET(1).EQ.0) THEN



!C.. Check if we're blocking the hamiltonian
!C      IF(THFBASIS.AND.TBLOCK) THEN
!C         WRITE(6,*) "THFBASIS set and NBLK=0.  ",
!C     &         "Cannot block diagonalize in HF Basis."
!C         STOP
!C      ENDIF
!C      CALL SYMGENEXCITS(FDET,NEL,NBASIS)
!C      CALL LeaveMemoryManager
!C      STOP


!C.. in order to calculate the H matrix, we need to work out all the determinants
!C.. beware with NPATHS - it calcs the list of dets even if we don't calc H
!C.. Could be big.
!C..Now we see how many determinants we need
!C      IF(nBasis.GT.170) THEN
!C..This fix is to stop floating overflow as taking the factorial of (nBasis.GT.170) crashes
!C  using the old FACTRL routine.
         NDET=1
         DNDET=1.0_dp
         DO I=0,NEL-1
            NDET=(NDET*(nBasis-I))/(I+1)
            DNDET=(DNDET*real(nBasis-I,dp))/real(I+1,dp)
         ENDDO

        IF (abs(real(NDET) - dndet) > 1.0e-6) THEN
!         WRITE(6,*) ' NUMBER OF DETERMINANTS : ' , DNDET
         NDET=-1
        ELSE
!         WRITE(6,*) ' NUMBER OF DETERMINANTS : ' , NDET
        ENDIF

!C      CALL TC(I_HMAX,I_P,NWHTAY)


!Check that the symmetry routines have set the symmetry up correctly...
        tSuccess=.true.
        tFoundOrbs(:)=.false.

        IF((.not.tHub).and.(.not.tUEG).and.TwoCycleSymGens) THEN
            do i=1,nSymLabels
!                WRITE(6,*) "NSymLabels: ",NSymLabels,i-1
                EndSymState=SymLabelCounts(1,i)+SymLabelCounts(2,i)-1
!                WRITE(6,*) "Number of states: ",SymLabelCounts(2,i)
                do j=SymLabelCounts(1,i),EndSymState

                    Beta=(2*SymLabelList(j))-1
                    Alpha=(2*SymLabelList(j))
                    SymAlpha=INT((G1(Alpha)%Sym%S),4)
                    SymBeta=INT((G1(Beta)%Sym%S),4)
!                    WRITE(6,*) "***",Alpha,Beta

                    IF(.not.tFoundOrbs(Beta)) THEN
                        tFoundOrbs(Beta)=.true.
                    ELSE
                        CALL Stop_All("SetupParameters","Orbital specified twice")
                    ENDIF
                    IF(.not.tFoundOrbs(Alpha)) THEN
                        tFoundOrbs(Alpha)=.true.
                    ELSE
                        CALL Stop_All("SetupParameters","Orbital specified twice")
                    ENDIF

                    IF(G1(Beta)%Ms.ne.-1) THEN
                        tSuccess=.false.
                    ELSEIF(G1(Alpha)%Ms.ne.1) THEN
                        tSuccess=.false.
                    ELSEIF((SymAlpha.ne.(i-1)).or.(SymBeta.ne.(i-1))) THEN
                        tSuccess=.false.
                    ENDIF
                enddo
            enddo
            do i=1,nBasis
                IF(.not.tFoundOrbs(i)) THEN
                    WRITE(6,*) "Orbital: ",i, " not found."
                    CALL Stop_All("SetupParameters","Orbital not found")
                ENDIF
            enddo
        ENDIF
        ! SpinOrbSymSetup currently sets up the symmetry arrays for use with
        ! symrandexcit2 excitation routines. These are not currently
        ! compatible with non-abelian symmetry groups, which CPMD jobs
        ! invariably used. To avoid this complication, this symmetry
        ! setup will not be used with CPMD, and thus these excitation
        ! generators won't work.
        IF(.not.tCPMD) THEN
            IF(.not.tSuccess) THEN
                WRITE(6,*) "************************************************"
                WRITE(6,*) "**                 WARNING!!!                 **"
                WRITE(6,*) "************************************************"
                WRITE(6,*) "Symmetry information not set up correctly in NECI initialisation"
                WRITE(6,*) "Will attempt to set up the symmetry again, but now in terms of spin orbitals"
                WRITE(6,*) "Old excitation generators will not work"
                WRITE(6,*) "I strongly suggest you check that the reference energy is correct."
                !CALL SpinOrbSymSetup() !.true.)
            ELSE
                WRITE(6,*) "Symmetry and spin of orbitals correctly set up for excitation generators."
                WRITE(6,*) "Simply transferring this into a spin orbital representation."
                !CALL SpinOrbSymSetup() !.false.)
            ENDIF

            if (tCSF) then
                call csf_sym_setup ()
            elseif (tPickVirtUniform) then
                call virt_uniform_sym_setup ()
            else
                ! Includes normal & HPHF
                call SpinOrbSymSetup ()
            endif

        ENDIF
! From now on, the orbitals are also contained in symlabellist2 and symlabelcounts2.
! These are stored using spin orbitals.

        ! If we are using CSFs, then we need to convert this into a csf
        if (tCSF) then
            ncsf = det_to_random_csf (FDET)
            write (6, '("Generated starting CSF: ")', advance='no')
            call write_det (6, FDET, .true.)
        endif

    End Subroutine DetInit

    function get_helement_compat (nI, nJ, IC, iLutI, iLutJ) result (hel)
        ! Get the matrix element of the hamiltonian. This assumes that we
        ! already know IC. We do not need to know iLutI, iLutJ (although
        ! they are helpful). This better fits the requirements of existing
        ! code than get_helement_normal.
        !
        ! In:  nI, nJ       - The determinants to consider
        !      iLutI, iLutJ - Bit representations of I,J (optional, helpful)
        !      IC           - The number of orbitals I,J differ by
        ! Ret: hel          - The desired matrix element.

        integer, intent(in) :: nI(nel), nJ(nel)
        integer(kind=n_int), intent(in), optional :: iLutI(0:NIfTot), iLutJ(0:NIfTot)
        integer, intent(in) :: IC
        HElement_t(dp) :: hel

        character(*), parameter :: this_routine = 'get_helement_compat'

        integer :: temp_ic

        if (tHPHFInts) &
            call stop_all (this_routine, "Should not be calling HPHF &
                          &integrals from here.")

        ! nobody actually uses Simons old CSF implementations..
!         if (tCSF) then
!             if (iscsf(nI) .or. iscsf(nJ)) then
!                 hel = CSFGetHelement (nI, nJ)
!                 return
!             endif
!         endif

        if (t_lattice_model) then
            temp_ic = ic
            hel = get_helement_lattice(nI, nJ, temp_ic)
            return
        end if

        if (tCSF) then
            if (iscsf(nI) .or. iscsf(nJ)) then
                hel = CSFGetHelement (nI, nJ)
                return
            endif
        endif

        if (tStoreAsExcitations) &
            call stop_all(this_routine, "tStoreExcitations not supported")

        if (present(iLutJ)) then
            hel = sltcnd_knowIC (nI, iLutI, iLutJ, IC)
        else
            hel = sltcnd_compat (nI, nJ, IC)
        endif

        ! Add in ECore if for a diagonal element
        if (IC == 0) then
            hel = hel + (ECore)
        else if (modk_offdiag) then
            hel = -abs(hel)
        end if

    end function

    function get_helement_normal (nI, nJ, iLutI, iLutJ, ICret) result(hel)

        ! Get the matrix element of the hamiltonian.
        !
        ! In:  nI, nJ       - The determinants to consider
        !      iLutI, iLutJ - Bit representations of I,J (optional, helpful)
        ! Out: ICret        - The number of orbitals I,J differ by
        ! Ret: hel          - The desired matrix element.

        integer, intent(in) :: nI(nel), nJ(nel)
        integer(kind=n_int), intent(in), optional :: iLutI(0:NIfTot), iLutJ(0:NIfTot)
        integer, intent(out), optional :: ICret
        HElement_t(dp) :: hel

        character(*), parameter :: this_routine = 'get_helement_normal'
        integer :: ex(2,2), IC
        integer(kind=n_int) :: ilut(0:NIfTot,2)

        if (tHPHFInts) &
            call stop_all (this_routine, "Should not be calling HPHF &
                          &integrals from here.")

        if (tCSF) then
            if (iscsf(nI) .or. iscsf(nJ)) then
                hel = CSFGetHelement (nI, nJ)
                return
            endif
        endif

        if (t_lattice_model) then
            if (present(ICret)) then
                ic = -1
                hel = get_helement_lattice(nI, nJ, ic)
                ICret = ic
            else
                hel = get_helement_lattice(nI,nJ)
            end if
            return
        end if

        if (tStoreAsExcitations .and. nI(1) == -1 .and. nJ(1) == -1) then
            ! TODO: how to express requirement for double?
            !if (IC /= 2) &
            !    call stop_all (this_routine, "tStoreAsExcitations in &
            !                  &get_helement requires IC=2 (doubles)")

            ex(1,:) = nJ(4:5)
            ex(2,:) = nJ(6:7)
            hel = sltcnd_2 (nI, ex, .false.)
        endif

        if (present(iLutJ)) then
            hel = sltcnd (nI, iLutI, iLutJ, IC)
        else
            call EncodeBitDet (nI, iLut(:,1))
            call EncodeBitdet (nJ, iLut(:,2))
            ! TODO: This is not an ideal place to end up...
            hel = sltcnd (nI, iLut(:,1), ilut(:,2), IC)
        endif

        ! Add in ECore for a diagonal element
        if (IC == 0) then
            hel = hel + (ECore)
        else if (modk_offdiag) then
            hel = -abs(hel)
        end if

        ! If requested, return IC
        if (present(ICret)) then
            ICret = IC
        endif

    end function get_helement_normal

    function get_helement_excit (NI, NJ, IC, ExcitMat, TParity) result(hel)

        ! Calculate the Hamiltonian Matrix Element for a given determinant (or
        ! csf), when we have the excitation matrix and parity of the
        ! excitation.
        !
        ! In:  nI, nJ       - The determinants to evaluate
        !      IC           - The number of orbitals I,J differ by
        !      ex           - The excitation matrix
        !      tParity      - The parity of the excitation
        ! Ret: hel          - The H matrix element

        integer, intent(in) :: nI(nel), nJ(nel), IC
        integer, intent(in) :: ExcitMat(2,ic)
        logical, intent(in) :: tParity
        HElement_t(dp) :: hel

        character(*), parameter :: this_routine = 'get_helement_excit'

        ! intermediately put the special call to the hubbard matrix elements
        ! here. Although I want to change that in the whole code to have
        ! procedure pointers similar to the excitation generator, which gets
        ! intialized to the correct function at the beginning of the
        ! excitations
        ! store all the lattice model matrix elements in one call.
        if (t_lattice_model) then
            hel = get_helement_lattice(nI, ic, ExcitMat, tParity)
            return
        end if

        ! If we are using CSFs, then call the csf routine.
        ! TODO: Passing through of ExcitMat to CSFGetHelement
        if (tCSF) then
            if (iscsf(NI) .or. iscsf(NJ)) then
                hel = CSFGetHelement (nI, nJ)
                return
            endif
        endif

        if (IC < 0) &
            call stop_all(this_routine, "get_helement_excit should only be &
                         &used if we know the number of excitations and the &
                         &excitation matrix")

        hel = sltcnd_excit_old(nI, IC, ExcitMat, tParity)

        if (IC == 0) then
            hel = hel + (ECore)
        else if (modk_offdiag) then
            hel = -abs(hel)
        end if

    end function get_helement_excit

    function get_helement_det_only (nI, nJ, iLutI, iLutJ, ic, ex, tParity, &
                                    HElGen) result (hel)

        ! Calculate the Hamiltonian Matrix Element for a determinant as above.
        ! This function assumes that we have got it correct for determinants
        ! (i.e. no error checking), and no conditionals. It also has extra
        ! arguments for compatibility with the function pointer methods.
        !
        ! In:  nI, nJ       - The determinants to evaluate
        !      iLutI, iLutJ - Bit representations (unused)
        !      ic           - The number of orbitals i,j differ by
        !      ex           - Excitation matrix
        !      tParity      - Parity of the excitation
        ! Ret: hel          - The H matrix element

        integer, intent(in) :: nI(nel), nJ(nel), ic, ex(2,ic)
        integer(kind=n_int), intent(in) :: iLutI(0:NIfTot), iLutJ(0:NIfTot)
        logical, intent(in) :: tParity
        HElement_t(dp) :: hel
        HElement_t(dp) , intent(in) :: HElGen    !Not used - here for compatibility with other interfaces.

        unused_var(ilutJ); unused_var(ilutI); unused_var(nJ); unused_var(hElgen);

        ! switch to lattice matrix element:
        if (t_lattice_model) then
            hel = get_helement_lattice(nI,ic, ex, tParity)
            return
        end if

        hel = sltcnd_excit_old(nI, IC, ex, tParity)

        if (IC == 0) then
            hel = hel + ECore
        else if (modk_offdiag) then
            hel = -abs(hel)
        end if
    end function

    HElement_t(dp) function GetH0Element4(nI,HFDet)
        ! Returns the matrix element of the unperturbed Hamiltonian,
        ! which is just the sum of the eigenvalues of the occupied
        ! orbitals and the core energy.
        ! HOWEVER, this routine is *SLOWER* than the GetH0Element3
        ! routine, and should only be used if you require this
        ! without reference to the fock eigenvalues.
        ! This is calculated by subtracting the required two electron terms
        ! from the diagonal matrix elements.
        integer, intent(in) :: nI(NEl),HFDet(NEl)
        HElement_t(dp) :: hel

        hel = SumFock(nI,HFDet)
        GetH0Element4 = hel + ECore

    end function GetH0Element4

    HElement_t(dp) function GetH0Element3(nI)
       ! Wrapper for GetH0Element.
       ! Returns the matrix element of the unperturbed Hamiltonian, which is
       ! just the sum of the eigenvalues of the occupied orbitals and the core
       ! energy.
       !  Note that GetH0Element{1,2} don't exist. The name is to be
       !  consistent with GetHElement3, i.e. offer the most abstraction possible.
       ! In:
       !    nI(nEl)  list of occupied spin orbitals in the determinant.
       integer nI(nEl)
       HElement_t(dp) hEl
       call GetH0Element(nI,nEl,Arr(1:nBasis,1:2),nBasis,ECore,hEl)
       GetH0Element3=hEl
    end function

    Subroutine DetCleanup()
    End Subroutine DetCleanup

    subroutine write_bit_rep(iUnit, iLut, lTerm)
       implicit none
       integer iUnit
       logical lTerm
       integer(n_int), intent(in) :: iLut(0:NIfTot)
       integer :: nI(nel), flags,i
       real(dp) :: sgn(lenof_sign)
       call extract_bit_rep(iLut,nI,sgn,flags)
       call write_det(iUnit,nI,.false.)
       write(iUnit,"(A)",advance='no') "("
       do i=1,lenof_sign
          write(iUnit, "(f16.7)", advance='no') sgn(i)
          if (i /= lenof_sign) write(iUnit, "(A)", advance='no') ","
       enddo
       write(iUnit,"(A,I5)", advance='no') ") ",flags
       if(lTerm) write(iUnit,*)
    end subroutine write_bit_rep

    subroutine get_lexicographic_dets (ilut_src, store, ilut_gen) !, det)

        integer(n_int), intent(in) :: ilut_src(0:NIfTot)
        type(lexicographic_store), intent(inout) :: store
        integer(n_int), intent(out), optional :: ilut_gen(0:NIfTot)
        !integer, intent(out), optional :: det(nel)

        integer :: i, nfound, orb, clro
        integer(n_int) :: ilut_tmp(0:NIfTot)

        ! If we haven't initialised the generator, do that now.
        if (.not. associated(store%dorder)) then

            ! Allocate dorder storage
            allocate(store%dorder(nel))
            store%dorder(1) = -1

            ! How many unpaired electrons are there
            store%nopen = count_open_orbs (ilut_src)
            store%nup = (store%nopen + LMS) / 2

            ! Obtain a list of unpaired orbitals
            nfound = 0
            !nelec = 0
            !allocate(store%open_indices(nopen))
            allocate(store%open_orbs(store%nopen))
            ilut_tmp = spatial_bit_det(ilut_src)
            do i = 1, nbasis-1, 2
                if (IsOcc(ilut_tmp, i)) then
                    if (IsOcc(ilut_tmp, i+1)) then
                    !    nelec = nelec + 2
                    else
                        nfound = nfound + 1
                    !    nelec = nelec + 1
                        store%open_orbs(nfound) = i
                    !    nhoce%open_indices(nfound) = nelec
                        if (nfound == store%nopen) exit
                    endif
                endif
            enddo
        endif

        ! Generate the next term in the sequence
        call get_lexicographic (store%dorder, store%nopen, store%nup)

        if (store%dorder(1) == -1) then
            deallocate(store%dorder)
            nullify(store%dorder)
            !deallocate(store%open_indices)
            deallocate(store%open_orbs)
            nullify(store%open_orbs)
            if (present(ilut_gen)) ilut_gen = 0
            !if (present(det)) det = 0
        else
            ! TODO: Test this with Ms /= 0.
            if (present(ilut_gen)) then
                ilut_gen = ilut_src
                do i = 1,store%nopen
                    if (store%dorder(i) == 0) then
                        orb = get_alpha(store%open_orbs(i))
                    else
                        orb = get_beta(store%open_orbs(i))
                    endif
                    set_orb(ilut_gen, orb)
                    clro=ab_pair(orb)
                    clr_orb(ilut_gen, clro)
                enddo
            endif

            !if (present(det)) ...
        endif

    end subroutine

END MODULE Determinants

      subroutine GetH0Element(nI,nEl,Arr,nBasis,ECore,hEl)
         !  Get a matrix element of the unperturbed Hamiltonian.  This is just
         !  the sum of the Hartree-Fock eigenvalues and the core energy.
         !  In:
         !     nI(nEl)  list of occupied spin orbitals in the determinant.
         !     nEl      # of electrons.
         !     Arr      array containing the eigenvalues of the spin-orbitals.
         !              (See System for how it's defined/stored).
         !     nBasis   # spin orbitals.
         !     ECore    Core energy.
         !  Out:
         !     hEl      <D_i|H_0|D_i>, the unperturbed Hamiltonian matrix element.
         use SystemData , only : TSTOREASEXCITATIONS
         use constants, only: dp
         implicit none
         integer nEl,nI(nEl),nBasis
         HElement_t(dp) hEl
         real(dp) Arr(nBasis,2),ECore
         integer i
         if(tStoreAsExcitations.and.nI(1).eq.-1) then
!The excitation storage starts with -1.  The next number is the excitation level,L .
!Next is the parity of the permutation required to lineup occupied->excited.  Then follows
!a list of the indexes of the L occupied orbitals within the HFDET, and then L virtual spinorbitals.
            hEl=0.0_dp
            do i=4,nI(2)+4-1
               hEl=hEl-(Arr(nI(i),2))
            enddo
            do i=i,i+nI(2)-1
               hEl=hEl+(Arr(nI(i),2))
            enddo
         else
            hEl=ECore
            do i=1,nEl
               hEl=hEl+(Arr(nI(i),2))
            enddo
         endif
!         call writedet(77,nI,nel,.false.)
!         write(77,*) "H0",hEl
!         call neci_flush(77)
      end subroutine

      subroutine DetFreezeBasis(GG)
        Use Determinants, only: FDet, nUHFDet, write_det, write_det_len
        use SystemData, only : nEl, nBasis, nBasisMax,BasisFN,G1,tFixLz
        use IntegralsData, only : nFrozen,nFrozenIn
        use sort_mod
        use sym_mod
        implicit none
        integer i,j
        INTEGER GG(*),Lz
        Type(BasisFn) s
        character(*), parameter :: this_routine = 'DetFreezeBasis'
!C.. Deal with FDET
!C.. GG(I) is the new position in G of the (old) orb I
         IF(FDET(1).NE.0) THEN
            J=0
            DO I=1,NEL
               FDET(I)=GG(FDET(I))
!C.. any orbitals which no longer exist, we move outside the basis
               IF(FDET(I).EQ.0) THEN
                  FDET(I)=nBasis+1
               ELSE
                  J=J+1
               ENDIF
            ENDDO
            call sort(fdet)
            IF(J.NE.NEL-NFROZEN-NFROZENIN) THEN
               WRITE(6,*) "Failed Freezing Det:"
               call write_det (6, FDET, .true.)
               call stop_all(this_routine, "After Freezing, FDET has wrong number of electrons")
            ENDIF
         ENDIF
         IF(nUHFDet(1).NE.0) THEN
            J=0
            DO I=1,NEL
               nUHFDET(I)=GG(nUHFDET(I))
!C.. any orbitals which no longer exist, we move outside the basis
               IF(nUHFDET(I).EQ.0) THEN
                  nUHFDET(I)=nBasis+1
               ELSE
                  J=J+1
               ENDIF
            ENDDO
            call sort (nUHFDet(1:nel))
            IF(J.NE.NEL-NFROZEN-NFROZENIN) THEN
               WRITE(6,*) "Failed Freezing Det:"
               call write_det (6, nUHFDET, .true.)
               call stop_all(this_routine, "After Freezing, UHFDET has wrong number of electrons")
            ENDIF
         ENDIF
         WRITE(6,"(A)",advance='no') " Post-Freeze Fermi det (D0):"
         call write_det_len (6, fDet, nel-nfrozen-nfrozenin, .true.)
         WRITE(6,"(A)",advance='no') " Symmetry: "
         Call GetSym(FDet,nEl-nFrozen-nFrozenIn,G1,nBasisMax,s)
         Call WriteSym(6,s%Sym,.true.)
         IF(tFixLz) THEN
             Call GetLz(FDet,nEl-nFrozen-nFrozenIn,Lz)
             WRITE(6,"(A,I5)") " Lz of Fermi det:",Lz
         ENDIF

      end subroutine


      LOGICAL FUNCTION ISUHFDET(NI,NEL)
         use SystemData , only : TUSEBRILLOUIN
         Use Determinants, only : NUHFDET
         IMPLICIT NONE
         INTEGER NEL,NI(NEL)
         INTEGER I
         IF(.NOT.TUSEBRILLOUIN) THEN
             ISUHFDET=.FALSE.
             RETURN
         ENDIF
         ISUHFDET=.TRUE.
         DO I=NEL,1,-1
            IF(NI(I).NE.NUHFDET(I)) THEN
                ISUHFDET=.FALSE.
                EXIT
            ENDIF
         ENDDO
!         ISUHFDET=.FALSE.
         RETURN
      END Function


! Generate the active space from a basis.
! The Active basis can be used to in PATHS calculations and (?) as a CASCI

! nActiveBasis(1:2) contains (First Active Basis Fn, Last Active Basis Fn)
! nDown is the number of orbital sets  below the Fermi level
! nUp is the number of orbital sets  above the Fermi level

      SUBROUTINE GenActiveBasis(ARR,nBasis,nEl,nActiveBasis, nDown,nUp)
         use SystemData, only: BasisFN
         use constants, only: dp
         IMPLICIT NONE
         INTEGER nEl,nActiveBasis(2),nBasis
         real(dp) ARR(nBasis)
         INTEGER I,nDown,nUp,nLeft
         I=nEl+1
         nLeft=1+nUp
         IF(nDown.NE.0.AND.nUp.NE.0) WRITE(6,*) "Including ",-nDown,",",nUp," extra degenerate sets in active space."
         DO WHILE (nLeft.GT.0.AND.I.LT.nBasis)
            DO WHILE (I.LT.nBasis.AND.ABS(ARR(I)-ARR(I-1)).LT.1.0e-5_dp)
               I=I+1
            ENDDO
            nLeft=nLeft-1
            IF(nLeft.EQ.nUp.AND.I.NE.nEl+1) WRITE(6,*) "Fermi determinant degenerate.  "
            IF(nLeft.ne.0) I=I+2
         ENDDO
         IF(I.EQ.nEl+1.and.nDown.eq.0) THEN
!We have no degeneracies at the Fermi Energy
            WRITE(6,*) "Fermi determinant non-degenerate.  "
            IF(nDown.eq.0) THEN
               WRITE(6,*) "Active space empty."
               nActiveBasis(1)=nEl+1
               nActiveBasis(2)=nEl
               RETURN
            ENDIF
         ENDIF
         nActiveBasis(2)=I-1
         I=nEl-1
         nLeft=nDown
         Do WHILE(nLeft.GT.0.AND.I.Gt.0)

            DO WHILE (I.GT.0.AND.ABS(ARR(I)-ARR(I+1)).LT.1.0e-5_dp)
               I=I-1
            ENDDO
            nLeft=nLeft-1
         ENDDO
         nActiveBasis(1)=I+1
         WRITE(6,*) "Active space:", nActiveBasis(1)," TO ",nActiveBasis(2)," (ordered labels)."
         WRITE(6,*) "Active space electrons:",nEl-nActiveBasis(1)+1
         RETURN
      END

      SUBROUTINE GENRANDOMDET(NEL,NBASIS,MCDET)
         use sort_mod
         use constants, only: dp
         IMPLICIT NONE
         INTEGER NEL,NBASIS,MCDET(NEL)
         INTEGER I,J,EL,SEED
         LOGICAL BR
         real(dp) RAN2
         SEED=-7
         DO I=1,NEL
            BR=.TRUE.
            DO WHILE (BR)
               BR=.FALSE.
               EL=INT(RAN2(SEED)*NBASIS+1)
               DO J=1,I-1
                  IF(MCDET(J).EQ.EL) BR=.TRUE.
               ENDDO
            ENDDO
            MCDET(I)=EL
         ENDDO
         call sort (mcDet)
         RETURN
      END

! Write determinant NI(NEL) to unit NUnit.  Set LTerm if to add a newline at end.  Also prints CSFs
      !SUBROUTINE WRITEDET(NUNIT,NI,NEL,LTERM)
      !   use legacy_data, only: CSF_NBSTART
      !   IMPLICIT NONE
      !   INTEGER NUNIT,NEL,NI(NEL),I
      !   LOGICAL LTERM
      !   INTEGER IEL
      !   CHARACTER*2 SUFF
      !   WRITE(NUNIT,"(A)",advance='no') "("
      !   DO I=1,NEL
      !      IEL=NI(I)
      !      IF(IEL.GE.CSF_NBSTART) THEN
      !         WRITE(NUNIT,"(I3)",advance='no'),(IEL-CSF_NBSTART)/4+1
      !         IEL=IAND(IEL-CSF_NBSTART,3)
      !         IF(IEL.EQ.0) THEN
      !            WRITE(NUNIT,"(A)",advance='no') "-B,"
      !         ELSEIF(IEL.EQ.1) THEN
      !            WRITE(NUNIT,"(A)",advance='no') "-A,"
      !         ELSEIF(IEL.EQ.2) THEN
      !            WRITE(NUNIT,"(A)",advance='no') "+B,"
      !         ELSE
      !            WRITE(NUNIT,"(A)",advance='no') "+A,"
      !         ENDIF
      !      ELSE
      !         WRITE(NUNIT,"(I5,A)",advance='no') IEL,","
      !      ENDIF
      !   ENDDO
      !   WRITE(NUNIT,"(A)",advance='no') ")"
      !   IF(LTERM) WRITE(NUNIT,*)
      !   RETURN
      !END
    subroutine writedet_oldcsf (nunit, nI, nel, lTerm)
        use systemdata, only: tCSF, tCSFOLD

        ! Write a human readable determinant to specified file unit. For use
        ! with old csf routines.
        ! Not easy to test, as both iscsf routines will return true sometimes.
        ! This is here for use if it becomes necessary (eg debugging)
        !
        ! In: nunit    - File unit
        !     nI (nel) - Determinant to print
        !     nel      - Number of electrons
        !     lTerm    - Do we write an end-of-line character

        use legacy_data, only: CSF_NBSTART
        implicit none
        integer, intent(in) :: nunit, nel, nI(nel)
        logical, intent(in) :: lTerm
        integer :: i, orb
        logical iscsf_old, bCSF

        ! Is this a csf? Note use of old (non-modularised) iscsf
        bCSF = (tCSF .or. tCSFOLD) .and. iscsf_old(nI, nel)

        ! Begin with an open bracket
        write(nunit,'("(")',advance='no')
        do i=1,nel
            orb = nI(i)
            if (bCSF .and. orb > CSF_NBSTART) then
                ! Output the component orbital
                orb = (orb - CSF_NBSTART) / 4 + 1
                write(nunit,'(i4)',advance='no') orb

                ! Output components of Yamanouchi symbol
                orb = nI(i) - csf_nbstart
                if (btest(orb, 1)) then
                    write(nunit,'("+")',advance='no')
                else
                    write(nunit,'("-")',advance='no')
                endif

                ! Output components of Ms
                if (btest(orb, 0)) then
                    write(nunit,'("A")',advance='no')
                else
                    write(nunit,'("B")',advance='no')
                endif
            else
                write(nunit,'(i6)',advance='no') orb
            endif
            if (i /= nel) write(nunit,'(",")',advance='no')
        enddo

        ! Close the written determinant off
        write(nunit,'(")")',advance='no')
        if (lTerm) write(nunit,*)

    end subroutine


! Calculate the one-electron part of the energy of a det
      FUNCTION CALCT(NI,NEL)
         use constants, only: dp
         USE SystemData, only : BasisFN
         USE OneEInts, only : GetTMatEl
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),I
         LOGICAL ISCSF_old
         real(dp) :: CALCT
         CALCT=0.0_dp
         IF(ISCSF_old(NI,NEL)) RETURN
         DO I=1,NEL
            CALCT=CALCT+GetTMATEl(NI(I),NI(I))
         ENDDO
         RETURN
      END

! Write bit-determinant NI to unit NUnit.  Set LTerm if to add a newline at end.  Also prints CSFs
      SUBROUTINE WriteBitDet(nUnit,iLutnI,lTerm)
         use SystemData, only : nEl
         use bit_reps, only: nIfTot
         use bit_reps, only: decode_bit_det
         use Determinants, only: write_det
         use constants, only: n_int
         implicit none
         integer nUnit,nI(nEl)
         integer(kind=n_int) :: iLutnI(0:nIfTot)
         logical lTerm
         call decode_bit_det (nI, iLutnI)
         call write_det (nUnit, nI, lTerm)
      END

      !Write out the determinant in bit notation
      SUBROUTINE WriteDetBit(nUnit,iLutnI,lTerm)
         use SystemData,only: nBasis
         use bit_reps, only: nIfTot
         use constants, only: n_int,bits_n_int
         implicit none
         integer, intent(in) :: nUnit
         integer(kind=n_int), intent(in) :: iLutnI(0:NIfTot)
         logical, intent(in) :: lTerm
         integer :: i

         do i=1,nBasis-1
             if(IsOcc(iLutnI,i)) then
                 write(nUnit,"(A1)",advance='no') "1"
             else
                 write(nUnit,"(A1)",advance='no') "0"
             endif
         enddo
         if(IsOcc(iLutnI,nBasis)) then
             if(lTerm) then
                 write(nUnit,"(A1)") "1"
             else
                 write(nUnit,"(A1)",advance='no') "1"
             endif
         else
             if(lTerm) then
                 write(nUnit,"(A1)") "0"
             else
                 write(nUnit,"(A1)",advance='no') "0"
             endif
         endif

      END SUBROUTINE WriteDetBit

! Write bit-determinant NI to unit NUnit.  Set LTerm if to add a newline at end.  Also prints CSFs
      SUBROUTINE WriteBitEx(nUnit,iLutRef,iLutnI,lTerm)
         use SystemData, only : nEl
         use bit_reps, only: NIfTot
         use constants, only: n_int
         implicit none
         integer nUnit
         integer(kind=n_int) :: iLutRef(0:nIfTot),iLutnI(0:nIfTot)
         integer Ex(2,nEl)
         logical lTerm
         logical tSign
         INTEGER iEl,I
         EX(1,1)=nEl  !Indicate the length of EX
         CALL GetBitExcitation(iLutRef,iLutnI,Ex,tSign)
         WRITE(NUNIT,"(A)",advance='no') "("
! First the excit from
         DO I=1,NEL
            IEL=EX(1,I)
            if(iEl.eq.0) EXIT
            WRITE(NUNIT,"(I5,A)",advance='no') IEL,","
         ENDDO
         IF(tSign) THEN
            WRITE(NUNIT,"(A)",advance='no') ")->-("
         ELSE
            WRITE(NUNIT,"(A)",advance='no') ")->+("
         ENDIF
!Now the excit to
         DO I=1,NEL
            IEL=EX(2,I)
            if(iEl.eq.0) EXIT
            WRITE(NUNIT,"(I5,A)",advance='no') IEL,","
         ENDDO
         WRITE(NUNIT,"(A)",advance='no') ")"
         IF(LTERM) WRITE(NUNIT,*)
         RETURN
      END

! This takes a the ground state FDet generated for the UEG and changes its total momentum
! According to input options
      subroutine ModifyMomentum(FDet)
        use SystemData, only : nEl,G1,k_momentum,nBasis,tUEG
        implicit none
        integer :: i,j ! Loop variables
        integer, intent(inout) :: FDet(NEl)
        integer :: k_total(3) ! Stores the total momentum of FDet
        integer :: delta_k(3) ! Stores the difference between the current FDet and the FDet we're aiming for
        integer, allocatable :: kPointToBasisFn(:,:,:,:) ! Look up table for kPoints to basis functions
        integer :: kmaxX,kminX,kmaxY,kminY,kmaxZ,kminZ,iSpinIndex ! Stores the limits of kPointToBasisFn
        integer :: det_sorted(NEl), e_store ! Storage for the sorting routine
        logical :: sorted ! As above
        integer :: wrapped_index
        integer :: k_new

        IF(.not.tUEG) call stop_all("ModifyMomentum", "Only works for UEG")

        ! Finds current momentum, and finds the difference between this and the target momentum
        ! most commonly will be zero to start with
        k_total(1)=0
        k_total(2)=0
        k_total(3)=0
        do j=1,NEl
            k_total(1)=k_total(1)+G1(FDet(j))%k(1)
            k_total(2)=k_total(2)+G1(FDet(j))%k(2)
            k_total(3)=k_total(3)+G1(FDet(j))%k(3)
        enddo
        delta_k=k_momentum-k_total

        if (delta_k(1).eq.0.and.delta_k(2).eq.0.and.delta_k(3).eq.0) write(6,*) "WARNING: specified momentum is ground state"

        ! Creates a look-up table for k-points (this was taken from symrandexcit2.F90)
        kmaxX=0
        kminX=0
        kmaxY=0
        kminY=0
        kminZ=0
        kmaxZ=0
        do i=1,nBasis
            IF(G1(i)%k(1).gt.kmaxX) kmaxX=G1(i)%k(1)
            IF(G1(i)%k(1).lt.kminX) kminX=G1(i)%k(1)
            IF(G1(i)%k(2).gt.kmaxY) kmaxY=G1(i)%k(2)
            IF(G1(i)%k(2).lt.kminY) kminY=G1(i)%k(2)
            IF(G1(i)%k(3).gt.kmaxZ) kmaxZ=G1(i)%k(3)
            IF(G1(i)%k(3).lt.kminZ) kminZ=G1(i)%k(3)
        enddo
        allocate(kPointToBasisFn(kminX:kmaxX,kminY:kmaxY,kminZ:kmaxZ,2))
        do i=1,nBasis
            iSpinIndex=(G1(i)%Ms+1)/2+1 ! iSpinIndex equals 1 for a beta spin (ms=-1), and 2 for an alpha spin (ms=1)
            kPointToBasisFn(G1(i)%k(1),G1(i)%k(2),G1(i)%k(3),iSpinIndex)=i
        enddo

        ! For each of the three dimensions, nudge electrons one at a time by one momentum unit until delta_k is reached
        det_sorted=FDet

        ! Bubble sort to order det_sorted in order of kx of the corresponding electron
        do
            sorted=.true.
            do i=1,NEl-1
                j=i+1
                if (G1(det_sorted(j))%k(1).gt.G1(det_sorted(i))%k(1)) then
                    sorted=.false.
                    e_store=det_sorted(i)
                    det_sorted(i)=det_sorted(j)
                    det_sorted(j)=e_store
                endif
            enddo
            if (sorted) exit
        enddo

        ! Nudge momenta one at a time
        if (delta_k(1).gt.0) then
            do i=1,delta_k(1)
                wrapped_index=mod(i,NEl)        ! Take the modulus to know which electron to nudge
                if (wrapped_index.eq.0) then    ! Deal with the i=NEl case
                    wrapped_index=NEl
                endif
                j=wrapped_index                 ! For convenience asign this to j
                k_new=G1(det_sorted(j))%k(1)+1  ! Find the new momentum of this electron
                if (k_new.gt.kmaxX) then        ! Check that this momentum isn't outside the cell
                    call stop_all("ModifyMomentum", "Electron moved outside of the cell limits")
                endif
                iSpinIndex=(G1(j)%Ms+1)/2+1     ! Spin of the new orbital is the same as the old
! Finds basis number for the new momentum
                det_sorted(j)=kPointToBasisFn(k_new,G1(det_sorted(j))%k(2),G1(det_sorted(j))%k(3),iSpinIndex)
            enddo
        else if (delta_k(1).lt.0) then ! For the negative case, i must run through negative numbers
            do i=-1,delta_k(1),-1
                wrapped_index=mod(i,NEl)
                if (wrapped_index.eq.0) then
                    wrapped_index=-NEl
                endif
                j=NEl+wrapped_index+1 ! Now this goes through the list backward (wrapped_index is negative)
                k_new=G1(det_sorted(j))%k(1)-1 ! Find the new momentum of this electron, this time in the opposite direction
                if (k_new.lt.kminX) then ! Check the limits of the cell again
                    call stop_all("ModifyMomentum", "Electron moved outside of the cell limits")
                endif
                iSpinIndex=(G1(j)%Ms+1)/2+1 ! Spin of the new orbital is the same as the old
! Finds basis number for the new momentum
                det_sorted(j)=kPointToBasisFn(k_new,G1(det_sorted(j))%k(2),G1(det_sorted(j))%k(3),iSpinIndex)
            enddo
        endif

        FDet=det_sorted

        !====ky treated as kx above
        do
            sorted=.true.
            do i=1,NEl-1
                j=i+1
                if (G1(det_sorted(j))%k(2).gt.G1(det_sorted(i))%k(2)) then
                    sorted=.false.
                    e_store=det_sorted(i)
                    det_sorted(i)=det_sorted(j)
                    det_sorted(j)=e_store
                endif
            enddo
            if (sorted) exit
        enddo

        ! Nudge momenta one at a time
        if (delta_k(2).gt.0) then
            do i=1,delta_k(2)
                wrapped_index=mod(i,NEl)        ! Take the modulus to know which electron to nudge
                if (wrapped_index.eq.0) then    ! Deal with the i=NEl case
                    wrapped_index=NEl
                endif
                j=wrapped_index                 ! For convenience asign this to j
                k_new=G1(det_sorted(j))%k(2)+1  ! Find the new momentum of this electron
                if (k_new.gt.kmaxY) then        ! Check that this momentum isn't outside the cell
                    call stop_all("ModifyMomentum", "Electron moved outside of the cell limits")
                endif
                iSpinIndex=(G1(j)%Ms+1)/2+1     ! Spin of the new orbital is the same as the old
! Finds basis number for the new momentum
                det_sorted(j)=kPointToBasisFn(G1(det_sorted(j))%k(1),k_new,G1(det_sorted(j))%k(3),iSpinIndex)
            enddo
        else if (delta_k(2).lt.0) then ! For the negative case, i must run through negative numbers
            do i=-1,delta_k(2),-1
                wrapped_index=mod(i,NEl)
                if (wrapped_index.eq.0) then
                    wrapped_index=-NEl
                endif
                j=NEl+wrapped_index+1 ! Now this goes through the list backward (wrapped_index is negative)
                k_new=G1(det_sorted(j))%k(2)-1 ! Find the new momentum of this electron, this time in the opposite direction
                if (k_new.lt.kminY) then ! Check the limits of the cell again
                    call stop_all("ModifyMomentum", "Electron moved outside of the cell limits")
                endif
                iSpinIndex=(G1(j)%Ms+1)/2+1 ! Spin of the new orbital is the same as the old
! Finds basis number for the new momentum
                det_sorted(j)=kPointToBasisFn(G1(det_sorted(j))%k(1),k_new,G1(det_sorted(j))%k(3),iSpinIndex)
            enddo
        endif

        FDet=det_sorted

        !====kz treated as kx and ky above
        do
            sorted=.true.
            do i=1,NEl-1
                j=i+1
                if (G1(det_sorted(j))%k(3).gt.G1(det_sorted(i))%k(3)) then
                    sorted=.false.
                    e_store=det_sorted(i)
                    det_sorted(i)=det_sorted(j)
                    det_sorted(j)=e_store
                endif
            enddo
            if (sorted) exit
        enddo

        ! Nudge momenta one at a time
        if (delta_k(3).gt.0) then
            do i=1,delta_k(3)
                wrapped_index=mod(i,NEl)        ! Take the modulus to know which electron to nudge
                if (wrapped_index.eq.0) then    ! Deal with the i=NEl case
                    wrapped_index=NEl
                endif
                j=wrapped_index                 ! For convenience asign this to j
                k_new=G1(det_sorted(j))%k(3)+1  ! Find the new momentum of this electron
                if (k_new.gt.kmaxZ) then        ! Check that this momentum isn't outside the cell
                    call stop_all("ModifyMomentum", "Electron moved outside of the cell limits")
                endif
                iSpinIndex=(G1(j)%Ms+1)/2+1     ! Spin of the new orbital is the same as the old
! Finds basis number for the new momentum
                det_sorted(j)=kPointToBasisFn(G1(det_sorted(j))%k(1),G1(det_sorted(j))%k(2),k_new,iSpinIndex)
            enddo
        else if (delta_k(3).lt.0) then ! For the negative case, i must run through negative numbers
            do i=-1,delta_k(3),-1
                wrapped_index=mod(i,NEl)
                if (wrapped_index.eq.0) then
                    wrapped_index=-NEl
                endif
                j=NEl+wrapped_index+1 ! Now this goes through the list backward (wrapped_index is negative)
                k_new=G1(det_sorted(j))%k(3)-1 ! Find the new momentum of this electron, this time in the opposite direction
                if (k_new.lt.kminZ) then ! Check the limits of the cell again
                    call stop_all("ModifyMomentum", "Electron moved outside of the cell limits")
                endif
                iSpinIndex=(G1(j)%Ms+1)/2+1 ! Spin of the new orbital is the same as the old
                ! Finds basis number for the new momentum
                det_sorted(j)=kPointToBasisFn(G1(det_sorted(j))%k(1),G1(det_sorted(j))%k(2),k_new,iSpinIndex)
            enddo
        endif

        FDet=det_sorted

        ! Bubble sort to order FDet back into increasing order by number
        do
            sorted=.true.
            do i=1,NEl-1
                j=i+1
                if (FDet(j).lt.FDet(i)) then
                    sorted=.false.
                    e_store=FDet(i)
                    FDet(i)=FDet(j)
                    FDet(j)=e_store
                endif
            enddo
            if (sorted) exit
        enddo

        write(6,*) "Total momentum set to", k_momentum

      end subroutine ModifyMomentum
