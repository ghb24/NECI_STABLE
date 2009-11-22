! This is a random excitation generator for use with csfs.
! Generate using a normalised and calculable probability.
module GenRandSymExcitCSF
    use Systemdata, only: nel, NIftot, tNoSymGenRandExcits, G1, LMS, nbasis
    use SymExcitDataMod
    use SymData, only: TwoCycleSymGens
    use csf, only: csf_orbital_mask, csf_test_bit, csf_apply_random_yama
    use csf, only: get_num_csfs, csf_apply_yama, csf_get_yamas
    use mt95, only: genrand_real2
    use GenRandSymExcitNUMod, only: ClassCountInd
    use DetBitOps, only: EncodeBitDet
    use Parallel
    implicit none
    
contains
    subroutine GenRandSymCSFExcit (nI, iLut, nJ, pSingle, IC, ExcitMat,&
                                   exFlag, pGen, ClassCountDoubleOcc2, &
                                   ClassCountSingleOcc2, ClassCountUnocc2,&
                                   tFilled)
        integer, intent(in)    :: nI(nel), iLut(0:NIfTot), exFlag
        integer, intent(out)   :: nJ(nel), IC, ExcitMat(2,2)
        integer, intent(inout) :: ClassCountDoubleOcc2(ScratchSize)
        integer, intent(inout) :: ClassCountSingleOcc2(ScratchSize)
        integer, intent(inout) :: ClassCountUnocc2(ScratchSize)
        logical, intent(inout) :: tFilled
        real*8,  intent(in)    :: pSingle
        real*8,  intent(out)   :: pGen
        character(*), parameter   :: this_routine = 'GenRandSymExcitCSF'
        integer :: Attempts, i, nopen
        logical :: bSingle

        ! Count the open shell electrons
        call CalcOpenOrbs (iLut, nopen)

        ! If the array is not already populated, perform an O[N] operation to
        ! find the number of occupied alpha/beta electrons, and number of
        ! occupied e- of each symmetry class and spin.
        if (.not. tFilled) then
            if ((.not.TwoCycleSymgens) .and. (.not.tNoSymGenRandExcits)) then
                write(6,'("GenRandSymExcitCSF can only be used for molecular&
                          & systems")')
                write(6,'("This is because of difficulties with other &
                          &symmetries setup.")')
                write(6,'("If you want to use these excitation generators &
                          &then add NOSYMGEN to the input to ignor symmetry &
                          &while generating excitations.")')
                call flush(6)
                call stop_all(this_routine,"GenRandsymExcitCSF can only be &
                                  &used for molecular systems using symmetry")
            endif

            call ConstructClassCounts(nI, nel-nopen, ClassCountDoubleOcc2, &
                                      ClassCountSingleOcc2, ClasscountUnocc2)

            tFilled = .true.
        endif

        select case (ExFlag)
            case (1)
                IC = 1
            case default
                call stop_all (this_routine, "Unsupported excitation mode")
        end select

        if (IC == 1) then
            call CreateSingleExcit (nI, nJ, ClassCountDoubleOcc2, &
                                ClassCountSingleOcc2, ClassCountUnocc2, &
                                iLut, ExcitMat, nopen, pSingle, pGen)
        endif
    end subroutine

    ! TODO: Do we need to include tParity?
    subroutine CreateSingleExcit (nI, nJ, CCDbl, CCSgl, CCUn, iLut, ExcitMat,&
                                  nopen, pSingle, pGen)
        integer, intent(in)  :: nI(nel), iLut(0:NIfTot), nopen
        integer, intent(out) :: nJ(nel), Excitmat(2,2)
        integer, intent(in)  :: CCSgl (ScratchSize) ! ClassCountSingleOcc2
        integer, intent(in)  :: CCDbl (ScratchSize) ! ClassCountDoubleOcc2
        integer, intent(in)  :: CCUn (ScratchSize)  ! ClassCountUnocc2
        real*8,  intent(in)  :: pSingle
        real*8,  intent(out) :: pGen
        character(*), parameter :: this_routine = 'CreateSingleExcit'
        integer :: elecsWNoExcits, elec, orb, orb2, spn, ind, norbs, symEx
        integer :: lnopen, ncsf, i, sym_ind, nexcit
        real*8 :: r, S
        logical :: bSingle

        ! Check that this condition is not necessary!!1
        !if (tNoSingsPossible .or. tNoSymGenRandExcits) then
        !    print*, 'nosingspossible', tNoSingsPossible
        !    print*, 'nosymgenrandexcits', tNoSymGenRandExcits
        !    call stop_all (this_routine, "This condition should never be &
        !                   &broken? See symrandexcit2.F90")
        !endif

        lnopen = nopen

        ! Loop over e- pairs to count e- with no excitations
        elecsWNoExcits = 0
        do i=2,ScratchSize,2
            ! Disallowed from beta e- in doubles
            if (CCDbl(i) /= 0) then
                elecsWNoExcits = elecsWNoExcits + CCDbl(i)
            endif

            ! Only from singles if more singles or vacant with same sym
            if ((CCSgl(i)/=0) .and. (CCSgl(i)<2) .and. (CCUn(i)==0)) then
                elecsWNoExcits = elecsWNoExcits + CCSgl(i)
            endif

            ! Alpha e- from doubles to singles or vacancies
            if ((CCDbl(i-1)/=0) .and. (CCSgl(i)==0) .and. (CCUn(i)==0)) then
                elecsWNoExcits = elecsWNoExcits + CCDbl(i-1)
            endif
        enddo

        ! 250 attempts to pick randomly
        do i=1,250
            ! Pick an electron at random, and extract its orbital
            call genrand_real2(r)
            elec = int(nel*r) + 1
            orb = iand(nI(elec), csf_orbital_mask)

            ! Obtain the symmetry index and Ms indicator (1=alpha)
            spn = 2
            if (G1(orb)%Ms == 1) spn = 1
            sym_ind = ClassCountInd(spn, int(G1(orb)%Sym%S,4), G1(orb)%Ml)

            ! Is this electron in a singly or doubly occupied orbital?
            ! Test if there are any allowed excitations
            orb2 = ieor((orb-1), 1)
            if (btest(iLut(orb2/32), mod(orb2,32))) then
                bSingle = .false.
                if (spn == 2) then
                    nexcit = 0
                else
                    nexcit = CCSgl(sym_ind+1) + CCUn(sym_ind+1)
                endif
            else
                if (spn == 1) call stop_all(this_routine, "Invalid Spin")
                bSingle = .true.
                nexcit = CCSgl(sym_ind) + (CCUn(sym_ind)-1)
            endif
            if (nexcit /= 0) exit
        enddo

        if (i > 250) then
            write(6,'("Cannot find single excitation after 250 attempts")')
            call stop_all(this_routine, "Cannot find single excitation after &
                                        &250 attempts")
        endif

        ! If we are exciting a singly occupied e-, then flip its spin.
        symEx = sym_ind
        if (bSingle) symEx = symEx - 1

        ! Choose an (allowed) unoccupied orbital to excite to. Draw orbitals
        ! from the desired symmetry and spin until we find one unoccupied.
        ! This is method 2 from symrandexcit2.
        norbs = OrbClassCount(symEx)
        do i=1,250
            call genrand_real2(r)
            orb2 = int(norbs*r)
            ind = SymLabelCounts2(1,symEx) + orb2
            orb2 = SymLabelList2(ind)

            ! Cannot excite a single to itself
            if (bSingle .and. (orb2 == orb+1)) cycle

            ! If target available, then select it. If exciting to a vacant
            ! orbital, then select the beta version.
            if (.not.(btest(iLut((orb2-1)/32), mod(orb2-1,32)))) then
                if (.not.btest(iLut((orb2-2)/32),mod(orb2-2,32))) then
                    orb2 = orb2-1
                    if (.not.bSingle) lnopen = lnopen + 2
                else if (bSingle) then
                    lnopen = lnopen - 2
                endif
                exit
            endif
        enddo

        if (i > 250) then
            write(6,'("Cannot find an unoccupied orbital for a single")')
            write(6,'("excitation after 250 attempts.")')
            write(6,'("Desired symmetry of unoccupied orbital =",i3)') &
                int(G1(orb)%Sym%S, 4)
            write(6,'("Num. orbitals (of correct spin) in symmetry =",i4)') &
                norbs
            write(6,'("Number of orbitals to legitimately pick =",i4)') nexcit
            call writedet(6,nI,nel,.true.)
            call stop_all(this_routine, "Cannot find an unoccupied orbital &
                         &for a single excitation after 250 attempts.")
        endif

        nJ = iand(nI, csf_orbital_mask)
        ! ExcitMat is the index of the orbital to excite from, and the actual
        ! orbital to excite to
        ExcitMat(1,1) = elec
        ExcitMat(2,1) = orb2

        ncsf = 1
        call csf_find_excit_det (ExcitMat, nJ, iLut, 1, nopen, lnopen, ncsf)

        ! TODO: test and remove these comments
        !call FindExcitDet(ExcitMat, nJ, 1, tParity)

        ! Turn this into a csf & apply a random Yamanouchi symbol
        !nJ = ior(nJ, csf_test_bit)
        !call csf_apply_random_yama (nJ, lnopen, real(LMS/2), ncsf)

        ! Generation probability
        pGen = pSingle / real(nexcit * (nel - elecsWNoExcits) * ncsf) 
    end subroutine

    ! Generate three arrays indicating the number of orbitals of each
    ! possible symmetry which are doubly-, singly- and un-occupied.
    subroutine ConstructClassCounts (nI, nclosed, CCDbl, CCSgl, CCUn)
        integer, intent(in) :: nI(nel), nclosed
        integer, intent(out) :: CCDbl (ScratchSize) ! ClassCountDoubleOcc2
        integer, intent(out) :: CCSgl (ScratchSize) ! ClassCountSingleOcc2
        integer, intent(out) :: CCUn (ScratchSize)  ! ClassCountUnocc2
        character(*), parameter :: this_routine = 'ConstructClassCounts'
        integer :: i, orb, ind

        ! nb. Unoccupied array is produced from overall orbital array minus
        !     the occupied electrons
        CCDbl = 0
        CCSgl = 0
        CCUn = OrbClasscount
        call writeDet(6,nI,nel,.true.)
        if (tNoSymGenRandExcits) then
            ! TODO: Implement ConstructClassCounts for tNoGenRandExcits
            call stop_all (this_routine, 'Unimplemented')
        else
            ! First loop over the closed shell electrons
            do i = 1,nclosed-1, 2
                ! Place e- into ClassCountDoubleOcc, and remove from Unocc.
                ! ind(beta) = ind(alpha) + 1 --> Can do both in one step.
                orb = iand (nI(i), csf_orbital_mask)
                ind = ClassCountInd(1, int(G1(orb)%Sym%S,4), G1(orb)%Ml)
                CCDbl(ind:ind+1) = CCDbl(ind:ind+1) + 1
                CCUn(ind:ind+1) = CCUn(ind:ind+1) - 1
            enddo

            ! Now loop over the open shell electrons
            do i = nclosed+1, nel
                orb = iand (nI(i), csf_orbital_mask)
                ind = ClassCountInd(2, int(G1(orb)%Sym%S,4), G1(orb)%Ml)
                CCSgl(ind) = CCSgl(ind) + 1
                CCUn(ind) = CCUn(ind) - 1
            enddo
        endif
    end subroutine

    ! Generate a determinant for the excitation specified in ExcitMat
    ! Note that this ASSUMES that you have got the allowed csf excitations
    ! correctly (a='alpha', b='beta', really just the occupation of spacial
    ! orbitals, always filling 'beta' first. Not really spin symmetry here)
    !:
    ! __SINGLES__:
    ! double a -> single a,   double a -> vacant b
    ! single b -> single a,   single b -> vacant b
    subroutine csf_find_excit_det (ExcitMat, nJ, iLut, IC, nopen, nopen_new,&
                                   ncsf, yama)
        integer, intent(in) :: IC, nopen, nopen_new, iLut(0:nIfTot)
        integer, intent(in), optional :: yama(ncsf, nopen_new)
        integer, intent(inout) :: ncsf, nJ(ncsf,nel), ExcitMat(2,2)
        integer :: i, pos, exbeta(IC), exalpha(IC), sralpha(IC), srbeta(IC)
        integer :: ins(4), nclosed
        character(*), parameter :: this_routine = "csf_find_excit_det"

        nJ(1,:) = iand(nJ(1,:), csf_orbital_mask)
        exbeta(1:IC) = ibclr(ExcitMat(2,1:IC)-1,0)+1
        exalpha(1:IC) = ibset(ExcitMat(2,1:IC)-1,0)+1
        srbeta(1:IC) = ibclr(nJ(1,ExcitMat(1,1:IC))-1,0)+1
        sralpha(1:IC) = ibset(nJ(1,ExcitMat(1,1:IC))-1,0)+1
        nclosed = nel - nopen
        if (IC == 1) then
            !>>>!print*, 'attempting single excitation'
            ! Are we exciting from a doubly occupied orbital?
            if (ExcitMat(1,1) <= nel-nopen) then
                !>>>! print*, 'exciting from a doubly occupied orbital'
                ! Are we exciting to a singly occupied, or a vacant orbital
                if (btest(iLut((exbeta(1)-1)/32),mod(exbeta(1)-1,32))) then
                    !>>>!print*, 'exciting to single', sralpha(1), exalpha(1)
                    ! Find the index of the beta e- in the spacial orbital
                    ! we are exciting to.
                    do i=nclosed+1,nel
                        if (nJ(1,i) == exbeta(1)) exit
                    enddo
                    if (i > nel) call stop_all (this_routine, &
                                                "Could not find orbital")

                    ! Remove this single from the list, and place the beta
                    ! e- from the original double into it.
                    nJ(1,i:nel-1) = nJ(1,i+1:nel)
                    ins(1) = srbeta(1)
                    call int_list_merge (nJ(1,nclosed+1:nel),ins(1:1),nopen,1)

                    ! Add the new double 
                    nJ(1,ExcitMat(1,1)-1:nclosed-2) = &
                            nJ(1,Excitmat(1,1)+1:nclosed)
                    ins(1) = exbeta(1)
                    ins(2) = exalpha(1)
                    call int_list_merge (nJ(1,1:nclosed),ins(1:2),nclosed-2,2)
                else
                    !>>>!print*, 'exciting to vacant', sralpha(1), exbeta(1)
                    ! Create two singles, both of them are betas.
                    ! One is the remaining e- from the original double.
                    ins(1) = min(srbeta(1),exbeta(1))
                    ins(2) = max(srbeta(1),exbeta(1))
                    nJ(1,ExcitMat(1,1)-1:nel-2) = nJ(1,Excitmat(1,1)+1:nel)
                    call int_list_merge (nJ(1,nclosed-1:nel), ins(1:2), &
                                         nopen, 2)
                endif
            else
                !>>>! print*,'exciting from a single'
                ! Are we exciting to a singly occupied, or a vacant orbital
                if (btest(iLut((exbeta(1)-1)/32),mod(exbeta(1)-1,32))) then
                    ! Test exciting a beta -> alpha
                    ! Remove both orbitals from singles
                    pos = nel
                    do i=nel,nclosed+1,-1
                        if (i == ExcitMat(1,1)) cycle
                        if (nJ(1,i) == ExcitMat(2,1)) cycle
                        nJ(1,pos) = nJ(1,i)
                        pos = pos - 1
                    enddo
                    ! Add to the list of doubles.
                    ins(1) = exbeta(1)
                    ins(2) = ExcitMat(2,1)
                    call int_list_merge (nJ(1,1:nclosed+2),ins(1:2),nclosed,2)
                else
                    ! Test exciting a beta -> beta
                    ! Remove orbital from singles, and insert a new one.
                    nJ(1,ExcitMat(1,1):nel-1) = nJ(1,Excitmat(1,1)+1:nel)
                    inS(1) = ExcitMat(2,1)
                    call int_list_merge (nJ(1,nclosed+1:nel), ins(1:1), &
                                         nopen-1, 1)                    
                endif
            endif
        else
            call stop_all(this_routine, "Not yet implemented")
        endif

        ! Make this into a csf that iscsf would recognise
        nJ(1,:) = ibset(nJ(1,:), csf_test_bit)

        ! If we have specified a yamanouchi symbol(s) apply it. Otherwise
        ! we must pick a random one.
        if (present(yama)) then
            forall (i=2:ncsf) nJ(i,:) = nJ(1,:)
            do i=1,ncsf
                call csf_apply_yama (nJ(i,:), yama(i,:))
            enddo
        else
            call csf_apply_random_yama (nJ, nopen_new, real(LMS/2,8), ncsf)
        endif
    end subroutine

    ! TODO: This currently ignores the possibility of just changing
    !       the yamanouchi symbol and not the configuration.
    ! TODO: Change from LMS --> STOT?
    subroutine csf_gen_excits (nI, iLut, nopen, bDouble, bSingle, CCDbl, &
                               CCSgl, CCUn, nexcit, nJ)
        integer, intent(in) :: nI(nel), ilut(0:NIfTot), nopen
        integer, intent(in) :: CCDbl(ScratchSize), CCSgl(ScratchSize)
        integer, intent(in) :: CCUn(ScratchSize)
        logical, intent(in) :: bDouble, bSingle
        integer, intent(out) :: nexcit
        integer, intent(out), dimension(:,:), allocatable, optional :: nJ
        character(*), parameter :: this_routine = 'csf_gen_excits'
        integer :: i, j, ind, sym_ind, spn, orb, orb2, numcsfs(-1:1), excit
        integer :: ierr, orb3, ExcitMat(2,2)
        integer, allocatable :: csf0 (:,:), csfp (:,:), csfm (:,:)
        real*8 :: S

        ! Calculate number of different Yamanouchi symbols given S
        ! and the possible values of nopen
        S = real(LMS) / 2
        numcsfs(0) = get_num_csfs (nopen, S)
        if (nopen<nel-1) numcsfs(1) = get_num_csfs (nopen+2, S)
        if (nopen>1) numcsfs(-1) = get_num_csfs (nopen-2, S)

        nexcit = 0
        if (bSingle) then
            ! Iterate over all the electrons. Select those with allowed
            ! transitions and sum the possible transitions
            do i=1,nel
                ! Obtain the orbital and its Ms/symmetry values
                orb = iand(nI(i), csf_orbital_mask)
                spn = (3 - G1(orb)%Ms) / 2 ! alpha=1, beta=2
                sym_ind = ClassCountInd(spn, int(G1(orb)%Sym%S,4), G1(orb)%Ml)

                ! Is it doubly or singly occupied
                orb2 = ieor((orb-1), 1)
                if (btest(iLut(orb2/32), mod(orb2,32))) then
                    ! Only allow transitions from doubly occupied alpha
                    if (spn == 1) then
                        nexcit = nexcit + (numcsfs(0)*CCSgl(sym_ind+1))
                        nexcit = nexcit + (numcsfs(1)*CCUn(sym_ind+1))
                    endif
                else
                    ! Only beta electrons allowed for singly occupied
                    if (spn == 1) call stop_all(this_routine, "Invalid spin")
                    nexcit = nexcit + numcsfs(-1)*(CCSgl(sym_ind)-1)
                    nexcit = nexcit + (numcsfs(0)*CCUn(sym_ind))
                endif
            enddo
        endif

        if (bDouble) then
            call stop_all(this_routine, 'Doubles not yet implemented')
        endif

        if (present(nJ)) then
            ! Allocate the required memory
            allocate(nJ(nexcit,nel), csf0(numcsfs(0),nopen), stat=ierr)
            if ((ierr == 0) .and. (nopen < nel-1)) &
                allocate(csfp(numcsfs(1), nopen+2), stat=ierr)
            if ((ierr == 0) .and. (nopen > 1)) & 
                allocate(csfm(numcsfs(-1), nopen-1), stat=ierr)
            if (ierr /= 0) call stop_all(this_routine, "Allocation failed")
            forall (i=1:nexcit) nJ(i,:) = nI

            ! Get all the required csfs
            call csf_get_yamas (nopen, S, csf0, numcsfs(0))
            if (nopen<nel-1) call csf_get_yamas (nopen+2, S, csfp, numcsfs(1))
            if (nopen>1) call csf_get_yamas (nopen-2, S, csfm, numcsfs(-1))

            ! Generate all the allowed singles
            excit = 1
            if (bSingle) then
                do i=1,nel
                    if (excit > nexcit) &
                        call stop_all(this_routine, "Generated too many csfs")

                    ! Obtain the orbital/symmetry to excite from
                    orb = iand(nI(i), csf_orbital_mask)
                    spn = (3 - G1(orb)%Ms) / 2 ! alpha=1, beta=2
                    sym_ind = ClassCountInd(spn, int(G1(orb)%Sym%S,4), &
                                            G1(orb)%Ml)
                    ExcitMat(1,1) = i
                                            
                    ! Is the source orbital doubly occupied?
                    orb2 = ieor((orb-1), 1) ! Spacial pair (zero based)
                    if (btest(iLut(orb2/32), mod(orb2,32))) then
                        ! Only promote alpha e- from doubly occupied orbitals
                        if (spn /= 1) cycle

                        ! Loop through all symmetry related orbitals
                        ind = SymLabelCounts2(1,sym_ind)
                        do j=1,OrbClassCount(sym_ind)
                            ! If the target orbital is filled, skip it
                            orb2 = SymLabelList2(ind+j-1)
                            if (btest(iLut((orb2-1)/32),mod(orb2-1,32))) cycle

                            ! Is this a vacant spacial orbital, or a single
                            orb3 = ieor((orb2-1),1) ! zero based
                            if (.not.btest(iLut(orb3/32),mod(orb3,32))) then
                                ! Excite into beta orbital of vacant pair
                                ExcitMat(2,1) = orb3+1
                                call csf_find_excit_det (ExcitMat, &
                                      nJ(excit:excit+numcsfs(1)-1,:), iLut, &
                                      1, nopen, nopen+2, numcsfs(1), csfp)
                                excit = excit + numcsfs(1)
                            else
                                ! Excite into alpha orbital of single
                                ExcitMat(2,1) = orb2
                                call csf_find_excit_det (ExcitMat, &
                                      nJ(excit:excit+numcsfs(0)-1,:), iLut, &
                                      1, nopen, nopen, numcsfs(0), csf0)
                                excit = excit + numcsfs(0)
                            endif
                        enddo
                    else ! Now consider excitations from singles.
                         ! Loop through all symmetry related orbitals
                         ! nb. spn == 2 (beta)
                         ind = SymLabelCounts2(1,sym_ind)
                         do j=1,OrbClassCount(sym_ind)
                             ! Cannot excite to self
                             orb2 = SymLabelList2(ind+j-1)
                             if (orb2 == orb) cycle

                             ! Is 'beta' orbital occipied?
                             if (btest(iLut((orb2-1)/32),mod(orb2-1,32))) then
                                 ! Check if 'alpha' is vacant (single->single)
                                 orb3 = ieor((orb2-1),1) ! zero based
                                 if (.not.btest(iLut(orb3/32),mod(orb3,32))) then
                                     ExcitMat(2,1) = orb3+1
                                     call csf_find_excit_det (ExcitMat, &
                                          nJ(excit:excit+numcsfs(-1)-1,:), iLut,&
                                          1, nopen, nopen-2, numcsfs(-1), csfm)
                                     excit = excit + numcsfs(-1)
                                 endif
                             else ! Exciting to vacant spacial pair (stay beta)
                                 ExcitMat(2,1) = orb2
                                 call csf_find_excit_det (ExcitMat, &
                                      nJ(excit:excit+numcsfs(0)-1,:), iLut, &
                                      1, nopen, nopen, numcsfs(0), csf0)
                                 excit = excit + numcsfs(0)
                             endif
                         enddo
                    endif
                enddo
            endif

            ! Clear up
            if (allocated(csf0)) deallocate (csf0)
            if (allocated(csfp)) deallocate (csfp)
            if (allocated(csfm)) deallocate (csfm)
        endif
    end subroutine

    subroutine TestCSF123 (nI)
        integer, intent(in) :: nI(nel)
        integer :: iLut(0:NIfTot), nopen
        integer :: CCDbl(ScratchSize), CCSgl(ScratchSize), CCUn(ScratchSize)
        integer :: ierr, nexcit, i
        integer, allocatable, dimension(:,:) :: nK
        character(*), parameter :: this_routine = 'TestGenRandSymCSFExcit'

        ! call TestGenRandSymCSFExcit (nI, 1000000, 1.0, 1, 10000)

        ! Generate bit representation, and count open shell electrons
        call EncodeBitDet (nI, iLut)
        call CalcOpenOrbs (iLut, nopen)

        print*, 'Starting determinant:'
        call writedet(6, nI, nel, .true.)

        ! Obtain the orbital symmetries for the following steps
        call ConstructClassCounts(nI, nel-nopen, CCDbl, CCSgl, CCUn)

        ! Enumerate all possible excitations
        call csf_gen_excits (nI, iLut, nopen, .false., .true., CCDbl, CCSgl,&
                             CCUn, nexcit, nK)
        print*, 'Excitations'
        do i=1,nexcit
            call writedet(6, nK(i,:), nel, .true.)
            call TestGenRandSymCSFExcit (nK(i,:), 1000000, 1.d0, 1,10000)
        enddo
        deallocate(nK)
    end subroutine

    ! A test routine for the CSF excitation generators. Initially generate
    ! (and count) all of the excited csfs. Then generate excititans randomly
    ! and histogram the generation probabilities.
    ! TODO: Doubles
    ! TODO: Only changing the Yamanouchi symbol
    subroutine TestGenRandSymCSFExcit (nI, iterations, pSingle, exFlag, &
                                       writeInterval)
        integer, intent(in) :: nI(nel), iterations, exFlag, writeInterval
        real*8,  intent(in) :: pSingle
        integer :: i, j, iLut(0:NIfTot), nJ(nel), ExcitMat(2,2), IC, nopen
        integer :: CCDbl(ScratchSize), CCSgl(ScratchSize), CCUn(ScratchSize)
        integer :: ierr, nexcit
        logical :: tFilled
        real*8  :: pGen, avContrib, avContribAll
        real*8,  allocatable, dimension(:,:) :: SinglesHist, AllSinglesHist
        integer, allocatable, dimension(:,:) :: nK
        character(*), parameter :: this_routine = 'TestGenRandSymCSFExcit'

        ! Generate bit representation, and count open shell electrons
        call EncodeBitDet (nI, iLut)
        call CalcOpenOrbs (iLut, nopen)

        ! Obtain the orbital symmetries for the following steps
        call ConstructClassCounts(nI, nel-nopen, CCDbl, CCSgl, CCUn)

        ! Enumerate all possible excitations
       ! call csf_gen_excits (nI, iLut, nopen, .false., .true., CCDbl, CCSgl,&
        !                     CCUn, nexcit, nK)
        !write(6,*), 'Excitations'
        !do i=1,nexcit
        !    call writedet(6, nK(i,:), nel, .true.)
        !enddo
        !deallocate(nK)
        ! If we don't want to generate them all, only count them.
        call csf_gen_excits (nI, iLut, nopen, .false., .true., CCDbl, CCSgl, &
                             CCUN, nexcit)

        ! Allocate memory for the histograms
        allocate (SinglesHist(nBasis,nBasis), &
                  AllSinglesHist(nBasis,nBasis), stat=ierr)
        if (ierr /= 0) call stop_all (this_routine,"Memory allocation failed")

        avContrib = 0
        avContribAll = 0
        SinglesHist = 0
        AllSinglesHist = 0
        tFilled = .true.
        open(9, file='AvContrib', status='unknown', position='append')
        do i=1,iterations
            ! Generate a random excitation
            call GenRandSymCSFExcit (nI, iLut, nJ, pSingle, IC, ExcitMat, &
                                     exFlag, pGen, CCDbl, CCSgl, CCUn,tFilled)

            ! Only average etc for an allowed transition
            if (nJ(1) /= 0) then
                avContrib = avContrib + 1/pGen

                if (IC == 1) then
                    SinglesHist(ExcitMat(1,1),ExcitMat(2,1)) = &
                        SinglesHist(ExcitMat(1,1),ExcitMat(2,1)) + (1/pGen)
                endif
            endif

            ! Take the average contribution over all processors, and write 
            ! out on node 0
            if (mod(i,writeInterval) == 0) then
                avContribAll = 0
#ifdef PARALLEL
                call MPI_Reduce (avContrib, avContribAll, 1, &
                                 MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
                                 MPI_COMM_WORLD, ierr)
#else
                AllAverageContrib = AverageContrib
#endif
                if (iProcIndex == 0) then
                    write(9,*) i, avContribAll/real(i*nexcit*nProcessors)
                endif
            endif
        enddo
        close(9)

#ifdef PARALLEL
        ! Sum the histograms over all processors
        call MPI_Reduce (SinglesHist, AllSinglesHist, nBasis**2, &
                         MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, &
                         ierr)
#else
        AllSinglesHist = SinglesHist
#endif

        ! Normalise the histograms and output in a readable form.
        ! These should tend to 0 or ncsf (an integer) for the excited csf.
        open (9,file="SinglesHist",status='unknown', position='append')
        do i=1,nbasis
            do j=1,nbasis
                if (AllSinglesHist(i,j) > 0) then
                    write(9,*)AllSinglesHist(i,j)/real(iterations*nProcessors)
                    pGen = AllSinglesHist(i,j)/real(iterations*nProcessors)
                    if (pGen > 10) then 
                        write (6,*) pGen,i,j,iand(ni(i),csf_orbital_mask)
                        call writedet(6,ni,nel,.true.)
                    endif                    
                endif
            enddo
        enddo
        close(9)

        ! Clean up
        deallocate (SinglesHist, AllSinglesHist)
    end subroutine
end module
