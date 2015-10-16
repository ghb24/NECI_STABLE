! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
module NatOrbsMod

    ! This file is primarily concerned with finding the one electron reduced
    ! density matrix, from the wavefunction constructed by a previous spawning
    ! calculation. Diagonalisation of this density matrix gives a set of
    ! eigenvectors which rotate the HF orbitals into the CI natural orbitals
    ! within the given excitation level. Once these eigenvectors have been
    ! obtained, the relevant routines from RotateOrbs are called to transform
    ! the integrals and produce a ROFCIDUMP file in the natural orbital basis.
    
    use Global_utilities
    use Parallel_neci
    use IntegralsData, only: UMAT
    use UMatCache, only: UMatInd
    use SystemData, only: NEl, nBasis, G1, ARR, BRR, lNoSymmetry, LMS, tStoreSpinOrbs, nOccAlpha, &
                          nOccBeta, tSeparateOccVirt, tRotateOccOnly, tRotateVirtOnly, &
                          tFindCINatOrbs, tUseMP2VarDenMat, nBasisMax, ALAT, iSpinSkip
    use bit_reps, only: NIfY, NIfTot
    use RotateOrbsData, only: SymLabelList2_rot, SymLabelCounts2_rot, SymLabelCounts2_rotTag, &
                              SymLabelListInv_rot, NoOrbs, SpatOrbs, FillOneRDM_time, &
                              FillMP2VDM_Time, DiagNatOrbMat_Time, OrderCoeff_Time, FillCoeff_Time, &
                              NoFrozenVirt, SymLabelList3_rot
    use sort_mod
    use bit_reps, only: decode_bit_det
    use MemoryManager, only: TagIntType
    use util_mod, only: get_free_unit
    use procedure_pointers, only: get_umat_el

    implicit none

    integer(TagIntType) :: NoSpinCyc, SymOrbs_rotTempTag
    real(dp), allocatable :: NatOrbMat(:,:), Evalues(:)
    integer, allocatable :: SymOrbs_rotTemp(:)
    integer(TagIntType) :: NatOrbMatTag, EvaluesTag

contains
    
    subroutine FindNatOrbs()

        ! Fed into this routine will be the wavefunction and its amplitudes
        ! within the given excitation level.    

        ! First need to set up the orbital labels and symmetries etc. This is
        ! done slightly differently for spin and spatial and whether or not we
        ! are truncating the virtual space when writing out the final ROFCIDUMP
        ! file.

        ! Allocate the matrix used to find the natural orbitals.

        integer :: ierr
        character(len=*), parameter :: t_r = 'FindNatOrbs'

        allocate(NatOrbMat(NoOrbs, NoOrbs),stat=ierr)
        call LogMemAlloc('NatOrbMat', NoOrbs**2,8,t_r, NatOrbMatTag,ierr)
        if (ierr /= 0) call stop_all(t_r,"Mem allocation for NatOrbMat failed.")
        NatOrbMat(:,:) = 0.0_dp

        allocate(Evalues(NoOrbs),stat=ierr)
        call LogMemAlloc('Evalues', NoOrbs,8,t_r,EvaluesTag,ierr)
        if (ierr /= 0) call stop_all(t_r,"Mem allocation for Evalues failed.")
        Evalues(:) = 0.0_dp

        ! First need to fill the relevant matrix for calculating the type of
        ! natural orbitals we want.
        if (tFindCINatOrbs) then

            ! For the CISD, CISDT etc natural orbitals, the relevant matrix is
            ! the one electron reduced density matrix from the previous
            ! spawning calculation (trucated at a certain excitation).
            call FillOneRDM()

        else if (tUseMP2VarDenMat) then

            ! For the MP2 natural orbitals, the natural orbital matrix is the
            ! MP2 variational density matrix.
            call FillMP2VDM()

        end if

        ! Then need to diagonalise this, maintaining the various symmetries
        ! (spin and spatial).
        call DiagNatOrbMat()

        ! We then need to put the resulting eigenvectors back into the ordering
        ! we want, and copy these over to CoeffT1.
        call OrderCoeffT1()
        
    end subroutine FindNatOrbs

    subroutine SetupNatOrbLabels()

        use MemoryManager, only: TagIntType

        integer :: x, i, j, ierr, NoOcc
        integer :: StartFill01, StartFill02, Symi, SymCurr, Prev, EndFill01, EndFill02
        character(len=*), parameter :: t_r = 'SetupNatOrbLabels'
        integer, allocatable :: LabVirtOrbs(:), LabOccOrbs(:), SymVirtOrbs(:), SymOccOrbs(:)
        integer(TagIntType) :: LabVirtOrbsTag, LabOccOrbsTag, SymVirtOrbsTag, SymOccOrbsTag
        integer :: lo, hi

        ! The earlier test should pick this up, if it crashes here, will want
        ! to put in an earlier test so that we don't get all the way to this
        ! stage.
        if ((LMS /= 0) .and. (.not.tStoreSpinOrbs)) then
            call stop_all("FindNatOrbs", "Open shell system, and UMAT is not being stored as spin orbitals.")
        end if

        ! We now need two slightly different sets of orbital labels for the
        ! case of spin orbitals and spatial orbitals. When using spin orbitals
        ! we want all the beta spin followed by all the alpha spin. Then we
        ! want two values for the number of occupied orbitals to allow for high
        ! spin cases. With spatial, it is equivalent to just keeping the beta
        ! spin.
        if (tStoreSpinOrbs) then
            NoSpinCyc = 2
        else
            NoSpinCyc = 1
        end if

        do x = 1, NoSpinCyc
            if (.not.tSeparateOccVirt) then
                NoOcc =  0
            else
                if (x == 1) then
                    if (tStoreSpinOrbs) then
                        NoOcc = nOccBeta
                    else
                        NoOcc = NEl/2
                    end if
                end if
                if (x == 2) NoOcc = nOccAlpha
            end if

            if (tSeparateOccVirt) then
                allocate(LabOccOrbs(NoOcc), stat=ierr)
                call LogMemAlloc('LabOccOrbs', (NoOcc), 4, t_r, LabOccOrbsTag, ierr)
                if (ierr /= 0) call stop_all(t_r, "Mem allocation for LabOccOrbs failed.")
                LabOccOrbs(:) = 0
                allocate(SymOccOrbs(NoOcc), stat=ierr)
                call LogMemAlloc('SymOccOrbs', (NoOcc), 4, t_r, SymOccOrbsTag, ierr)
                if (ierr /= 0) call stop_all(t_r, "Mem allocation for SymOccOrbs failed.")
                SymOccOrbs(:) = 0
            end if

            allocate(LabVirtOrbs(SpatOrbs-NoOcc), stat=ierr)
            call LogMemAlloc('LabVirtOrbs', (SpatOrbs-NoOcc), 4, t_r, LabVirtOrbsTag, ierr)
            if (ierr /= 0) call stop_all(t_r, "Mem allocation for LabVirtOrbs failed.")
            LabVirtOrbs(:) = 0
            allocate(SymVirtOrbs(SpatOrbs-NoOcc), stat=ierr)
            call LogMemAlloc('SymVirtOrbs', (SpatOrbs-NoOcc), 4, t_r, SymVirtOrbsTag, ierr)
            if (ierr /= 0) call stop_all(t_r, "Mem allocation for SymVirtOrbs failed.")
            SymVirtOrbs(:) = 0

            ! First fill SymLabelList2_rot.

            ! Brr has the orbital numbers in order of energy...
            ! i.e Brr(2) = the orbital index with the second lowest energy.

            ! This picks out the NoOcc lowest energy orbitals from BRR as these
            ! will be the occupied. These are then ordered according to
            ! symmetry, and the same done to the virtual.
            do i = 1, NoOcc
                if (x == 1) then
                    if (tStoreSpinOrbs) then
                        LabOccOrbs(i) = BRR(2*i)-1
                        SymOccOrbs(i) = int(G1(LabOccOrbs(i))%sym%S,4)
                    else
                        LabOccOrbs(i) = BRR(2*i)/2
                        SymOccOrbs(i) = int(G1(LabOccOrbs(i)*2)%sym%S,4)
                    end if
                else if (x == 2) then
                    LabOccOrbs(i) = BRR(2*i)
                    SymOccOrbs(i) = int(G1(LabOccOrbs(i))%sym%S,4)
                end if
            end do

            ! Sorts LabOrbs according to the order of SymOccOrbs (i.e. in
            ! terms of symmetry). 
            if (tSeparateOccVirt) call sort(SymOccOrbs, LabOccOrbs)

            ! Same for the virtual.
            do i = 1, SpatOrbs-NoOcc
                if (x == 1) then
                    if (tStoreSpinOrbs) then
                        if (tSeparateOccVirt) then
                            LabVirtOrbs(i) = BRR((2*i)+(NoOcc*2))-1
                            SymVirtOrbs(i) = int(G1(LabVirtOrbs(i))%sym%S,4)
                        else
                            LabVirtOrbs(i) = BRR((2*i))-1
                            SymVirtOrbs(i) = int(G1(LabVirtOrbs(i))%sym%S,4)
                        end if
                    else
                        if (tSeparateOccVirt) then
                            LabVirtOrbs(i) = BRR((2*i)+(NoOcc*2))/2
                            SymVirtOrbs(i) = int(G1(LabVirtOrbs(i)*2)%sym%S,4)
                        else
                            LabVirtOrbs(i) = BRR((2*i))/2
                            SymVirtOrbs(i) = int(G1(LabVirtOrbs(i)*2)%sym%S,4)
                        end if
                    end if
                else if (x == 2) then
                    if (tSeparateOccVirt) then
                        LabVirtOrbs(i) = BRR((2*i)+(NoOcc*2))
                        SymVirtOrbs(i) = int(G1(LabVirtOrbs(i))%sym%S,4)
                    else
                        LabVirtOrbs(i) = BRR((2*i))
                        SymVirtOrbs(i) = int(G1(LabVirtOrbs(i))%sym%S,4)
                    end if
                end if
            end do

            ! Sorts LabOrbs according to the order of SymOccOrbs (i.e. in
            ! terms of symmetry). 
            call sort(SymVirtOrbs, LabVirtOrbs)

            ! SymLabelList2_rot is then filled with the symmetry ordered
            ! occupied then virtual arrays for each spin.
            if (x == 1) then
                StartFill01 = 1
                StartFill02 = NoOcc+1
                EndFill01 = NoOcc
                EndFill02 = SpatOrbs
            else if (x == 2) then
                StartFill01 = SpatOrbs+1
                StartFill02 = SpatOrbs+NoOcc+1
                EndFill01 = SpatOrbs+NoOcc
                EndFill02 = NoOrbs
            end if

            j = 0
            do i = StartFill01,EndFill01
                j = j+1
                SymLabelList2_rot(i) = LabOccOrbs(j)
            end do
            j = 0
            do i = StartFill02, EndFill02
                j = j+1
                SymLabelList2_rot(i) = LabVirtOrbs(j)
            end do

            ! Second fill SymLabelCounts2_rot.
            ! - the first 8 places of SymLabelCounts2_rot(1,:) and
            !     SymLabelCounts2_rot(2,:) refer to the occupied orbitals 
            ! - and the second 8 to the virtuals.

            if (lNoSymmetry) then
                ! If we are ignoring symmetry, all orbitals essentially have
                ! symmetry 0.
                if (x == 1) then
                    SymLabelCounts2_rot(1,1) = 1
                    SymLabelCounts2_rot(1,9) = NoOcc+1
                    SymLabelCounts2_rot(2,1) = NoOcc
                    SymLabelCounts2_rot(2,9) = SpatOrbs-NoOcc
                else if (x == 2) then
                    SymLabelCounts2_rot(1,17) = 1
                    SymLabelCounts2_rot(1,25) = NoOcc+1
                    SymLabelCounts2_rot(2,17) = NoOcc
                    SymLabelCounts2_rot(2,25) = SpatOrbs-NoOcc
                end if
                
            else 
                ! Otherwise we run through the occupied orbitals, counting the
                ! number with each symmetry and noting where in
                ! SymLabelList2_rot each symmetry block starts.
                if (x == 1) then
                    StartFill01 = 1
                    StartFill02 = 9
                    Prev = 0
                else if (x == 2) then
                    StartFill01 = 17
                    StartFill02 = 25
                    Prev = SpatOrbs
                end if
                SymCurr = 0
                SymLabelCounts2_rot(1, StartFill01) = 1+Prev
                do i = 1, NoOcc
                    if (tStoreSpinOrbs) then
                        Symi = int(G1(SymLabelList2_rot(i+Prev))%sym%S,4)
                    else
                        Symi = int(G1((SymLabelList2_rot(i+Prev)*2))%sym%S,4)
                    end if
                    SymLabelCounts2_rot(2,(Symi+StartFill01)) = SymLabelCounts2_rot(2,(Symi+StartFill01))+1
                    if (Symi /= SymCurr) then
                        SymLabelCounts2_rot(1,(Symi+StartFill01)) = i+Prev
                        SymCurr= Symi
                    end if
                end do
                ! The same is then done for the virtuals.
                SymCurr = 0
                SymLabelCounts2_rot(1, StartFill02) = NoOcc+1+Prev
                do i = NoOcc+1, SpatOrbs
                    if (tStoreSpinOrbs) then
                        Symi = int(G1(SymLabelList2_rot(i+Prev))%sym%S,4)
                    else
                        Symi = int(G1((SymLabelList2_rot(i+Prev)*2))%sym%S,4)
                    end if
                    SymLabelCounts2_rot(2,(Symi+StartFill02)) = SymLabelCounts2_rot(2,(Symi+StartFill02))+1
                    if (Symi /= SymCurr) then
                        SymLabelCounts2_rot(1,(Symi+StartFill02)) = i+Prev
                        SymCurr = Symi
                    end if
                end do
            end if
     
            ! Go through each symmetry group, making sure the orbital pairs
            ! are ordered lowest to highest.
            if (x == 1) then
                do i = 1, 16
                    if (SymLabelCounts2_rot(2,i) /= 0) then
                        lo = SymLabelCounts2_rot(1, i)
                        hi = lo + SymLabelCounts2_rot(2, i) - 1
                        call sort(SymLabelList2_rot (lo:hi))
                    end if
                end do
            else if (x == 2) then
                do i = 17, 32
                    if (SymLabelCounts2_rot(2,i) /= 0) then
                        lo = SymLabelCounts2_rot(1, i)
                        hi = lo + SymLabelCounts2_rot(2, i) - 1
                        call sort(SymLabelList2_rot (lo:hi))
                    end if
                end do
            end if

            ! Deallocate the arrays just used in this routine.
            if (tSeparateOccVirt) then
                deallocate(SymOccOrbs)
                call LogMemDealloc(t_r, SymOccOrbsTag)

                deallocate(LabOccOrbs)
                call LogMemDealloc(t_r,LabOccOrbsTag)
            end if

            deallocate(SymVirtOrbs)
            call LogMemDealloc(t_r, SymVirtOrbsTag)

            deallocate(LabVirtOrbs)
            call LogMemDealloc(t_r,LabVirtOrbsTag)

        end do

        do i = 1, NoOrbs
            SymLabelListInv_rot(SymLabelList2_rot(i)) = i
        end do

        do i = 1, NoOrbs
            SymLabelList3_rot(i) = SymLabelList2_rot(i)
        end do

        if (.not.tSeparateOccVirt) then
            ! Basically we treat all the orbitals as virtuals and set NoOcc
            ! to zero in each routine. 
            tRotateVirtOnly = .true.
        end if

    end subroutine SetupNatOrbLabels

    subroutine FillOneRDM()

        ! Det is the number of determinants in FCIDets.
        ! FCIDets contains the list of all determinants in the system in bit
        ! string representation, FCIDets(0:NIfTot,1:Det)
        ! ICILevel is the max excitation level of the calculation - as in
        ! EXCITE ICILevel.
        ! FCIDetIndex(1:NEl) contains the index of FCIDets where each
        ! excitation level starts.
        ! As in FCIDetIndex(1) = 2 always I think - Excitation level 1 starts
        ! at the second determinant (after HF).
        ! Pretty sure FCIDetIndex always goes from 1:NEl even from truncated
        ! excite calculations.

        ! The elements of AllHistogram correspond to the rows of FCIDets - i.e
        ! to each determinant in the system. AllHistogram contains the final
        ! (normalised) amplitude of the determinant - with sign.

        use DetCalcData, only: Det,FCIDets,FCIDetIndex,ICILevel
        use DetBitOps, only: FindBitExcitLevel
        use bit_reps, only: decode_bit_det
        use hist_data, only: AllHistogram

        integer :: excit, i, j
        integer :: Starti, Endi, Startj, Endj, ExcitLevel, Ex(2,1)
        integer :: Orbi, Orbj, nJ(NEl), Orbk, k, nI(NEl), MaxExcit
        integer :: Spins
        logical :: tSign
        real(dp) :: SignDet

        ! Density matrix = D_pq = < Psi | a_q+ a_p | Psi > 
        !                = sum_ij [ c_i* c_j < D_i | a_q+ a_p | D_j > ]
        ! Where a_p is the annihilation, and a_q+ the creation operators.
        ! In other words, < D_i | a_q+ a_p | D_j > will only be non-zero if
        ! D_i and D_j are connected by an annihilation in p and a creation in q.

        ! Want to run through all determinants D_i in the final wavefunction.
        ! For each, find all determinants, D_j that are connected to D_i by a
        ! single excitation - i.e. those which differ by just one orbital. Only
        ! need to consider those of the same excitation level or one above or
        ! one below. Find the differing orbitals - these will be p and q. Sum
        ! in the occupations of D_i and D_j (c_i x c_j) to the matrix element
        ! D_pq. Take, for instance, p always =< q.

        ! Will get the orbitals in the original labelling system - convert it
        ! to this system.

        ! Get a list of the wavefunction with amplitudes in order of excitation.

        ! Depending on the type of reduced density matrix want to:
        ! Run through the determinants with excitation level one less, the same
        ! and one more.
        
        FillOneRDM_Time%timer_name = 'FillOneRDM'
        call set_timer(FillOneRDM_Time, 30)

        write(6,*) '*** The weight of the HF determinant is : ', AllHistogram(1,1)

        write(6,*) 'Beginning to fill the one-electron reduced density matrix.'

        if (ICILevel == 0) then
            MaxExcit = NEl
        else
            MaxExcit = ICILevel
        end if

        ! Run through all determinants D_i, in the final wavefunction, Psi. 
        ! If this is done by excitation block, we then don't have to check 
        ! the excitation level of the determinant each time.
        do excit = 0, MaxExcit         

            ! The HF only involves 'occupied' orbitals - these are not required
            ! if only rotating virt.
            if (tRotateVirtOnly .and. tSeparateOccVirt .and. (excit == 0)) cycle      

            ! This next bit is a bit messy because there is no row in
            ! FCIDetIndex for the HF - there is probably a tidier way to
            ! achieve the same thing, but it does the trick for now.
            if (excit == 0) then ! i is the HF det.
                Starti = 1
                Endi = 1
                Startj = 1
                Endj = min((FCIDetIndex(2)-1), Det) ! If i is the HF det, just run over singly excited j. 
            else if (excit == MaxExcit) then
                Starti = FCIDetIndex(excit)
                Endi = Det
                Startj = FCIDetIndex(excit-1)
                Endj = Det
            else
                Starti = FCIDetIndex(excit)
                Endi = FCIDetIndex(excit+1)-1
                if (excit == 1) then
                    Startj = 1
                    if (NEl < 3) then
                        Endj = Det
                    else
                        Endj = FCIDetIndex(3)-1
                    end if
                else if (excit == (MaxExcit-1)) then
                    Startj = FCIDetIndex(excit-1)
                    Endj = Det
                else
                    Startj = FCIDetIndex(excit-1)
                    Endj = FCIDetIndex(excit+2)-1
                end if
            end if

            ! Then run through the determinants in that excitation level.
            do i = Starti, Endi

                ! Run through all determinants D_j, with the potential to be
                ! connected to i by a single excitation, i.e from one
                ! excitation lower to one excitation higher.
                do j = Startj, i
                    if ((i > Det) .or. (j > Det)) then
                        call stop_all('FillOneRDM', 'Running through i or j larger than the number of determinants.')
                    end if

                    ExcitLevel = FindBitExcitLevel(FCIDets(:,i), FCIDets(:,j),2)

                    ! Need to find the excitation level between D_i and D_j. 
                    ! If this is 1 - go on to add their contributions to the
                    ! OneRDM.
                    if (ExcitLevel == 1) then
                        Ex(:,:) = 0
                        Ex(1,1) = ExcitLevel

                        call GetBitExcitation(FCIDets(:,i), FCIDets(:,j), Ex, tSign)
                        ! Gives the orbitals involved in the excitation Ex(1,1)
                        ! in i -> Ex(2,1) in j (in spin orbitals).

                        if (tStoreSpinOrbs) then
                            ! OneRDM will be in spin orbitals - simply add the
                            ! orbitals involved.
                            Orbi = SymLabelListInv_rot(Ex(1,1))
                            Orbj = SymLabelListInv_rot(Ex(2,1))
                            Spins = 1
                        else
                            Orbi = SymLabelListInv_rot(ceiling(real(Ex(1,1),dp)/2.0_dp))
                            Orbj = SymLabelListInv_rot(ceiling(real(Ex(2,1),dp)/2.0_dp))
                            Spins = 2
                        end if

                        if (tSign) then
                            SignDet = -1.0_dp
                        else
                            SignDet = 1.0_dp
                        end if

                        NatOrbMat(Orbi,Orbj) = NatOrbMat(Orbi,Orbj) + (SignDet*AllHistogram(1,i)*AllHistogram(1,j))
                        NatOrbMat(Orbj,Orbi) = NatOrbMat(Orbj,Orbi) + (SignDet*AllHistogram(1,i)*AllHistogram(1,j))

                        if ((AllHistogram(1,i)*AllHistogram(1,j) /= 0.0_dp) .and. &
                         (int(G1(SymLabelList2_rot(Orbi)*Spins)%sym%S,4) /= &
                         int(G1(SymLabelList2_rot(Orbj)*Spins)%sym%S,4))) then
                            write(6,*) 'ERROR in symmetries'
                            write(6,*) 'Ex,', Ex(1,1), Ex(2,1)
                            write(6,*) ceiling(real(Ex(1,1)/2.0_dp,dp)), ceiling(real(Ex(2,1)/2.0_dp,dp))
                            write(6,*) 'Orbi,', Orbi, 'Orbj,', Orbj
                            write(6,*) 'Sym(Orbi)', int(G1(SymLabelList2_rot(Orbi)*Spins)%sym%S,4),'Sym(Orbj)', &
                                int(G1(SymLabelList2_rot(Orbj)*Spins)%sym%S,4)
                            call decode_bit_det(nI, FCIDets(0:NIfTot,i))
                            write(6,*) 'i', nI
                            call decode_bit_det (nJ, FCIDets(0:NIfTot,j))
                            write(6,*) 'j', nJ
                            write(6,*) 'AllHistogram(1,i)', AllHistogram(1,i)
                            write(6,*) 'AllHistogram(1,j)', AllHistogram(1,j)
                            call neci_flush(6)
                            call stop_all('FillOneRDM','Non-zero element between different symmetries.')
                        end if

                    else if (ExcitLevel == 0) then
                        call Decode_Bit_Det(nJ,FCIDets(0:NIfTot,j))
                        do k = 1, NEl
                            if (tStoreSpinOrbs) then
                                Orbk = SymLabelListInv_rot(nJ(k))
                            else
                                Orbk = SymLabelListInv_rot(ceiling(real(nJ(k),dp)/2.0_dp))
                            end if
                            NatOrbMat(Orbk,Orbk) = NatOrbMat(Orbk,Orbk)+(AllHistogram(1,j)**2)
                            ! 0.5 x because this will be added twice since we
                            ! are not currently restricting i<k or anything.
                        end do
                    end if
                        
                end do

            end do
        end do

        write(6,*) 'Filled OneRDM'

        call halt_timer(FillOneRDM_Time)

    end subroutine FillOneRDM

    subroutine FillMP2VDM()

        ! In this routine, the natural orbital matrix is calculated from the
        ! MP2 variational density matrix.

        ! MP2VDM = D2_ab = sum_ijc [ t_ij^ac ( 2 t_ij^bc - t_ji^bc ) ]
        ! Where:  t_ij^ac = - < ab | ij > / ( E_a - E_i + E_b - Ej )
        ! Ref: J. Chem. Phys. 131, 034113 (2009) - note: in Eqn 1, the
        ! cb indices are the wrong way round (should be bc).

        use SystemData, only: tUEG
        use constants, only: dp

        integer :: a, b, c, i, j, a2, b2, c2, i2, j2, x, y, z, w
        integer :: Startab, Endab, NoOcc, NoOccC, Startc, Endc, Starti, Endi, Startj, Endj
        real(dp) :: MP2VDMSum
        character(len=*), parameter :: t_r = 'FillMP2VDM'
        HElement_t(dp) :: HEl01,HEl02

#ifdef __CMPLX
         call stop_all('FillMP2VDM', 'Natural Orbitals not implemented for complex orbitals.')
#endif
        ! Calculating the MP2VDM (D2_ab) matrix whose eigenvectors become the
        ! transformation matrix. This goes in the natural orbital matrix of
        ! this module. The eigenvalues are the occupation numbers of the new
        ! orbitals. These should decrease exponentially so that when we remove
        ! the orbitals with small occupation numbers we should have little
        ! effect on the energy.

        ! For the MP2VDM, we always only rotate the virtual orbitals -
        ! the denomonator term of the above expression would be 0 if a and b
        ! were occupied. The orbital labels are ordered occupied then virtual
        ! if spatial orbitals are being used, otherwise they go occupied beta,
        ! virtual beta, occupied alpha, virtual alpha. This is so the alpha and
        ! beta spins can be diagonalised separately and we can keep track of
        ! which is which when the evectors are reordered and maintain spin
        ! symmetry.

        write(6,*) 'Filling MP2VDM nat orb matrix'
        call neci_flush(6)
        
        FillMP2VDM_Time%timer_name = 'FillMP2VDM'
        call set_timer(FillMP2VDM_Time,30)

        do x = 1, NoSpinCyc
            if (x == 1) then
                if (tStoreSpinOrbs) then
                    NoOcc = nOccBeta
                else
                    NoOcc = NEl/2
                end if
                Startab = NoOcc+1
                Endab = SpatOrbs
            else if (x == 2) then
                NoOcc = nOccAlpha
                Startab = SpatOrbs+NoOcc+1
                Endab = NoOrbs
            end if

            ! a and b must be of the same spin to mix, so only need to run over
            ! both beta then both alpha.
            do a2 = Startab,Endab
                a = SymLabelList2_rot(a2)
                do b2 = Startab,a2

                    b = SymLabelList2_rot(b2)

                    MP2VDMSum = 0.0_dp

                    ! When a and b beta, run over both alpha and beta virtual
                    ! for c, then both alpha and beta virtual for both i and j
                    ! etc. 

                    do y = 1, NoSpinCyc
                        if (y == 1) then
                            if (tStoreSpinOrbs) then
                                NoOccC = nOccBeta
                            else
                                NoOccC = NEl/2
                            end if
                            Startc = NoOccC+1
                            Endc = SpatOrbs
                        else if (y == 2) then
                            NoOccC = nOccAlpha
                            Startc = SpatOrbs+NoOccC+1
                            Endc = NoOrbs
                        end if


                        do c2 = Startc,Endc
                            c = SymLabelList2_rot(c2)

                            do z = 1, NoSpinCyc
                                if (z == 1) then
                                    Starti = 1
                                    if (tStoreSpinOrbs) then
                                        Endi = nOccBeta
                                    else
                                        Endi = NEl/2
                                    end if
                                else if (z == 2) then
                                    Starti = 1+SpatOrbs
                                    Endi = SpatOrbs+nOccAlpha
                                end if

                                do i2 = Starti,Endi
                                    i = SymLabelList2_rot(i2)

                                    do w = 1, NoSpinCyc
                                        if (w == 1) then
                                            Startj = 1
                                            if (tStoreSpinOrbs) then
                                                Endj = nOccBeta
                                            else
                                                Endj = NEl/2
                                            end if
                                        else if (w == 2) then
                                            Startj = 1+SpatOrbs
                                            Endj = SpatOrbs+nOccAlpha
                                        end if


                                        do j2 = Startj, Endj
                                            j = SymLabelList2_rot(j2)

                                            if (tUEG) then
                                                HEl01 = get_umat_el(a, c, i, j)
                                                HEl02 = get_umat_el(b, c, i, j)
                                                MP2VDMSum = MP2VDMSum + &
                                                    &(( (real(HEl01,dp)) * (2.0_dp*(real(HEl02,dp))) )/&
                                                    &( (ARR(2*i,2)+ARR(2*j,2)-ARR(2*a,2)-ARR(2*c,2)) &
                                                    &* (ARR(2*i,2)+ARR(2*j,2)-ARR(2*b,2)-ARR(2*c,2)) ) )

                                                HEl02 = get_umat_el(c, b, i, j)
                                                MP2VDMSum = MP2VDMSum - &
                                                    &(( (real(HEl01,dp)) * (real(HEl02,dp)) )/&
                                                    &( (ARR(2*i,2)+ARR(2*j,2)-ARR(2*a,2)-ARR(2*c,2)) * &
                                                    &(ARR(2*i,2)+ARR(2*j,2)-ARR(2*c,2)-ARR(2*b,2)) ) )

                                            else if (tStoreSpinOrbs) then
                                                if ((ARR(i,2)+ARR(j,2)-ARR(a,2)-ARR(c,2)) == 0.0_dp) then
                                                    if ((real(UMAT(UMatInd(a,c,i,j,0,0)),dp)) /= 0.0_dp) then
                                                        write(6,*) i, j, a, c, real(UMAT(UMatInd(a, c, i, j, 0, 0)), dp)
                                                        call stop_all(t_r, "Dividing a non-zero by zero.")
                                                    end if
                                                end if
                                                MP2VDMSum = MP2VDMSum + &
                                                   (((real(UMAT(UMatInd(a, c, i, j, 0, 0)),dp)) & 
                                                   * (2.0_dp*(real(UMAT(UMatInd(b, c, i, j, 0, 0)),dp))))/&
                                                   ( (ARR(i,2)+ARR(j,2)-ARR(a,2)-ARR(c,2)) &
                                                   * (ARR(i,2)+ARR(j,2)-ARR(b,2)-ARR(c,2)) ) )
                                                MP2VDMSum=MP2VDMSum-&
                                                    (( (real(UMAT(UMatInd(a, c, i, j, 0, 0)),dp)) &
                                                    * (real(UMAT(UMatInd(c, b, i, j, 0, 0)),dp)) )/ &
                                                    ( (ARR(i,2)+ARR(j,2)-ARR(a,2)-ARR(c,2)) &
                                                    * (ARR(i,2)+ARR(j,2)-ARR(c,2)-ARR(b,2)) ) )
                                            else
                                                MP2VDMSum = MP2VDMSum + &
                                                    (( (real(UMAT(UMatInd(a, c, i, j, 0, 0)),dp)) &
                                                    * (2.0_dp*(real(UMAT(UMatInd(b, c, i, j, 0, 0)),dp))) )/&
                                                    ((ARR(2*i,2)+ARR(2*j,2)-ARR(2*a,2)-ARR(2*c,2)) &
                                                    * (ARR(2*i,2)+ARR(2*j,2)-ARR(2*b,2)-ARR(2*c,2))))
                                               MP2VDMSum = MP2VDMSum - &
                                                    (( (real(UMAT(UMatInd(a, c, i, j, 0, 0)),dp)) &
                                                    * (real(UMAT(UMatInd(c, b, i, j, 0, 0)),dp)) )/&
                                                    ( (ARR(2*i,2)+ARR(2*j,2)-ARR(2*a,2)-ARR(2*c,2)) &
                                                    * (ARR(2*i,2)+ARR(2*j,2)-ARR(2*c,2)-ARR(2*b,2))))
                                            end if

                                        end do
                                    end do
                                end do
                            end do
                        end do
                    end do
                    NatOrbMat(a2,b2) = MP2VDMSum
                    NatOrbMat(b2,a2) = MP2VDMSum
                end do
            end do
        end do

        write(6,*) 'Finished filling MP2VDM'

        call halt_timer(FillMP2VDM_Time)

    end subroutine FillMP2VDM

    subroutine DiagNatOrbMat()

        ! The diagonalisation routine reorders the orbitals in such a way that
        ! the corresponding orbital labels are lost. In order to keep the spin
        ! and spatial symmetries, each symmetry must be fed into the
        ! diagonalisation routine separately. The best way to do this is to
        ! order the orbitals so that all the alpha orbitals follow all the beta
        ! orbitals, with the occupied orbitals first, in terms of symmetry, and
        ! the virtual second, also ordered by symmetry. This gives us
        ! flexibility w.r.t rotating only the occupied or only virtual and
        ! looking at high spin states.

        use MemoryManager, only: TagIntType

        real(dp) :: SumTrace, SumDiagTrace
        real(dp), allocatable :: WORK2(:), EvaluesSym(:), NOMSym(:,:)
        integer :: ierr, i, j, x, z, Sym, LWORK2
        integer :: SymStartInd, NoSymBlock, PrevSym, StartOccVirt, EndOccVirt, Prev, NoOcc
        integer(TagIntType) :: EvaluesSymTag, NOMSymTag, WORK2Tag
        character(len=*), parameter :: t_r = 'DiagNatOrbMat'
 
        DiagNatOrbMat_Time%timer_name = 'DiagNatOrb'
        call set_timer(DiagNatOrbMat_Time, 30)

        do x = 1, NoSpinCyc
            if (tSeparateOccVirt) then
                if (x == 1) then
                    if (tStoreSpinOrbs) then
                        NoOcc = nOccBeta
                    else
                        NoOcc = NEl/2
                    end if
                    Prev = 0
                else if (x == 2) then
                    NoOcc = nOccAlpha
                    Prev = SpatOrbs
                end if
            else
                NoOcc = 0
            end if
            if (tRotateVirtOnly) then
                do i = 1, NoOcc
                    do j = 1, SpatOrbs
                        NatOrbMat(i+Prev,j+Prev) = 0.0_dp
                        NatOrbMat(j+Prev,i+Prev) = 0.0_dp
                        if (i == j) NatOrbMat(i+Prev, j+Prev) = 1.0_dp
                    end do
                    Evalues(i+Prev) = 1.0_dp
                end do
            else if (tRotateOccOnly) then
                do i = NoOcc+1, SpatOrbs
                    do j = 1, SpatOrbs
                        NatOrbMat(i+Prev, j+Prev) = 0.0_dp
                        NatOrbMat(j+Prev, i+Prev) = 0.0_dp
                        if (i == j) NatOrbMat(i+Prev, j+Prev) = 1.0_dp
                    end do
                    Evalues(i+Prev) = 1.0_dp
                end do
            else if (tSeparateOccVirt) then
                do i = 1, NoOcc
                    do j = NoOcc+1, SpatOrbs
                        NatOrbMat(i+Prev, j+Prev) = 0.0_dp
                        NatOrbMat(j+Prev, i+Prev) = 0.0_dp
                    end do
                end do
            end if
        end do

        ! Test that we're not breaking symmetry.
        do i = 1, NoOrbs
            do j = 1, NoOrbs
                if (tStoreSpinOrbs) then
                    if ((int(G1(SymLabelList2_rot(i))%sym%S,4) /= int(G1(SymLabelList2_rot(j))%sym%S,4))) then
                        if (abs(NatOrbMat(i,j)) >= 1.0e-15_dp) then
                            write(6,'(6A8,A20)') 'i', 'j', 'Label i', 'Label j', 'Sym i', 'Sym j', 'Matrix value'
                            write(6,'(6I3,F40.20)') i, j, SymLabelList2_rot(i), SymLabelList2_rot(j), &
                                int(G1(SymLabelList2_rot(i))%sym%S,4), &
                                int(G1(SymLabelList2_rot(j))%sym%S,4), NatOrbMat(i,j)
                            if (tUseMP2VarDenMat) then
                                write(6,*) '**WARNING** - There is a non-zero NatOrbMat value between " &
                                 & //"orbitals of different symmetry.'
                                write(6,*) 'These elements will be ignored, and the symmetry maintained " &
                                 & //"in the final transformation matrix.'
                            else
                                call stop_all(t_r, 'Non-zero NatOrbMat value between different symmetries.')
                            end if
                        end if
                        NatOrbMat(i,j) = 0.0_dp
                    end if
                else
                    if ((int(G1(SymLabelList2_rot(i)*2)%sym%S,4) /= int(G1(SymLabelList2_rot(j)*2)%sym%S,4))) then
                        if (abs(NatOrbMat(i,j)) >= 1.0e-15_dp) then
                            write(6,'(6A8,A20)') 'i', 'j', 'Label i', 'Label j', 'Sym i', 'Sym j', 'Matrix value'
                            write(6,'(6I3,F40.20)') i, j, SymLabelList2_rot(i), SymLabelList2_rot(j), &
                             int(G1(SymLabelList2_rot(i)*2)%sym%S,4), int(G1(SymLabelList2_rot(j)*2)%sym%S, 4), NatOrbMat(i,j)
                            if (tUseMP2VarDenMat) then
                                write(6,*) '**WARNING** - There is a non-zero NatOrbMat value between orbitals of " &
                                & //"different symmetry.'
                                write(6,*) 'These elements will be ignored, and the symmetry maintained in the " &
                                & //"final transformation matrix.'
                            else
                                call stop_all(t_r, 'Non-zero NatOrbMat value between different symmetries.')
                            end if
                        end if
                        NatOrbMat(i,j) = 0.0_dp
                    end if
                end if
            end do
        end do

        SumTrace = 0.0_dp
        do i = 1, NoOrbs
            SumTrace= SumTrace+NatOrbMat(i,i)
        end do

        write(6,*) 'Calculating eigenvectors and eigenvalues of NatOrbMat'
        call neci_flush(6)

        ! If we are using spin orbitals, need to feed in the alpha and beta
        ! spins separately. Otherwise these jumble up and the final ordering
        ! is uncorrect. There should be no non-zero values between these, but
        ! can put a check in for this.

        do x = 1, NoSpinCyc

            ! If we want to maintain the symmetry, we cannot have all the
            ! orbitals jumbled up when the diagonaliser reorders the
            ! eigenvectors. Must instead feed each symmetry block in
            ! separately. This means that although the transformed orbitals
            ! are jumbled within the symmetry blocks, the symmetry labels are
            ! all that are relevant and these are unaffected.
            StartOccVirt = 1
            EndOccVirt = 2
            if (tRotateVirtOnly) StartOccVirt = 2
            if (tRotateOccOnly) EndOccVirt = 1

            do z = StartOccVirt, EndOccVirt
                if (x == 1) then
                    if (z == 1) PrevSym = 1
                    if (z == 2) PrevSym = 9
                else if (x == 2) then
                    if (z == 1) PrevSym = 17
                    if (z == 2) PrevSym = 25
                end if

                Sym = 0
                LWORK2 = -1

                do while (Sym <= 7)

                    NoSymBlock = SymLabelCounts2_rot(2, Sym+PrevSym)

                    ! This is one less than the index that the symmetry starts,
                    ! so that when we run through i = 1,...,  we can start at
                    ! SymStartInd+i.
                    SymStartInd = SymLabelCounts2_rot(1, Sym+PrevSym)-1

                    if (NoSymBlock > 1) then

                        allocate(NOMSym(NoSymBlock, NoSymBlock), stat=ierr)
                        call LogMemAlloc('NOMSym', NoSymBlock**2, 8, t_r, NOMSymTag, ierr)
                        if (ierr /= 0) call stop_all(t_r, "Problem allocating NOMSym.")
                        allocate(EvaluesSym(NoSymBlock), stat=ierr)
                        call LogMemAlloc('EvaluesSym', NoSymBlock, 8, t_r, EvaluesSymTag, ierr)
                        if (ierr /= 0) call stop_all(t_r, "Problem allocating EvaluesSym.")

                        LWORK2 = 3*NoSymBlock + 1
                        allocate(WORK2(LWORK2), stat=ierr)
                        call LogMemAlloc('WORK2', LWORK2, 8, t_r, WORK2Tag, ierr)
                        if (ierr /= 0) call stop_all(t_r, "Problem allocating WORK2.")

                        do j = 1, NoSymBlock
                            do i = 1, NoSymBlock
                                NOMSym(i,j) = NatOrbMat(SymStartInd+i, SymStartInd+j)
                            end do
                        end do

                        write(6,*) '*****'
                        write(6,*) 'Symmetry ', Sym, 'with x ',x,' and z ',z,' has ', NoSymBlock,' orbitals.'
                        write(6,*) 'The NatOrbMat for this symmetry block is '
                        do i = 1, NoSymBlock
                            do j = 1, NoSymBlock
                                write(6,'(F20.10)',advance='no') NOMSym(j,i)
                            end do
                            write(6,*) ''
                        end do

                        ! NOMSym goes in as the original NOMSym, comes out as
                        ! the eigenvectors (Coefficients). EvaluesSym comes out
                        ! as the eigenvalues in ascending order.
                        call dsyev('V','L', NoSymBlock, NOMSym, NoSymBlock, EvaluesSym, WORK2, LWORK2, ierr)

                        write(6,*) 'After diagonalization, the e-vectors (diagonal elements) of this matrix are,'
                        do i = 1, NoSymBlock
                            write(6,'(F20.10)',advance='no') EvaluesSym(i)
                        end do
                        write(6,*) ''
                        write(6,*) 'These go from orbital,', SymStartInd+1,' to ', SymStartInd+NoSymBlock

                        do i = 1, NoSymBlock
                            Evalues(SymStartInd+i) = EvaluesSym(i)
                        end do

                        ! CAREFUL if eigenvalues are put in ascending order,
                        ! this may not be correct, with the labelling system.
                        ! Maybe better to just take coefficients and transform 
                        ! TMAT2DRot in transform2elints. A check that comes out
                        ! as diagonal is a check of this routine anyway.

                        write(6,*) 'The eigenvectors (coefficients) for symmetry block ', Sym
                        do i = 1, NoSymBlock
                            do j = 1, NoSymBlock
                                write(6,'(F20.10)',advance='no') NOMSym(j,i)
                            end do
                            write(6,*) ''
                        end do
                 
                        ! Directly fill the coefficient matrix with the
                        ! eigenvectors from the diagonalization.
                        do j = 1, NoSymBlock
                            do i = 1, NoSymBlock
                                NatOrbMat(SymStartInd+i, SymStartInd+j) = NOMSym(i,j)
                            end do
                        end do

                        deallocate(WORK2)
                        call LogMemDealloc(t_r,WORK2Tag)

                        deallocate(NOMSym)
                        call LogMemDealloc(t_r, NOMSymTag)

                        deallocate(EvaluesSym)
                        call LogMemDealloc(t_r,EvaluesSymTag)

                    else if (NoSymBlock == 1) then
                        ! The eigenvalue is the lone value, while the
                        ! eigenvector is 1.

                        Evalues(SymStartInd+1) = NatOrbMat(SymStartInd+1, SymStartInd+1)
                        NatOrbMat(SymStartInd+1, SymStartInd+1) =1.0_dp
                        write(6,*) '*****'
                        write(6,*) 'Symmetry ', Sym,' has only one orbital.'
                        write(6,*) 'Copying diagonal element,', SymStartInd+1,'to NatOrbMat'
                    end if

                    Sym= Sym+1
                end do
            end do
        end do

        write(6,*) 'Matrix diagonalised'
        call neci_flush(6)

        SumDiagTrace= 0.0_dp
        do i = 1, NoOrbs
            SumDiagTrace= SumDiagTrace+Evalues(i)
        end do
        if ((abs(SumDiagTrace-SumTrace)) > 10.0_dp) then
            write(6,*) 'Sum of diagonal NatOrbMat elements : ', SumTrace
            write(6,*) 'Sum of eigenvalues : ', SumDiagTrace
            write(6,*) 'WARNING, The trace of the 1RDM matrix before diagonalisation is not equal to that after.'
        end if

        call halt_timer(DiagNatOrbMat_Time)

    end subroutine DiagNatOrbMat

    subroutine OrderCoeffT1()

        use RotateOrbsData, only: SymLabelList3_rot
        use LoggingData, only: tTruncRODump

        integer :: x, i, ierr, StartSort, EndSort, NoOcc
        character(len=*), parameter :: t_r = 'OrderCoeffT1'
        
        ! Here, if symmetry is kept, we are going to have to reorder the
        ! eigenvectors according to the size of the eigenvalues, while taking
        ! the orbital labels (and therefore symmetries) with them. This will be
        ! put back into MP2VDM from MP2VDMTemp.

        ! Want to reorder the eigenvalues from largest to smallest, taking the
        ! eigenvectors with them and the symmetry as well. If using spin
        ! orbitals, do this for the alpha spin and then the beta.
 
        OrderCoeff_Time%timer_name = 'OrderCoeff'
        call set_timer(OrderCoeff_Time, 30)

        if (tTruncRODump) then
            ! If we are truncating, the orbitals stay in this order, so we want
            ! to take their symmetries with them.
            allocate(SymOrbs_rotTemp(NoOrbs), stat=ierr)
            call LogMemAlloc('SymOrbs_rotTemp', NoOrbs, 4, t_r, SymOrbs_rotTempTag, ierr)
            SymOrbs_rotTemp(:) = 0

            if (tStoreSpinOrbs) then
                do i = 1, NoOrbs
                    SymOrbs_rotTemp(i) = int(G1(SymLabelList2_rot(i))%sym%S, 4)
                end do
            else 
                do i = 1, NoOrbs
                    SymOrbs_rotTemp(i) = int(G1(SymLabelList2_rot(i)*2)%sym%S, 4)
                end do
            end if

            do x = 1, NoSpinCyc

                if (x == 1) then
                    if (tSeparateOccVirt) then
                        if (tStoreSpinOrbs) then
                            NoOcc = nOccBeta
                        else
                            NoOcc = NEl/2
                        end if
                    else
                        NoOcc = 0
                    end if
                    StartSort = 1
                    EndSort = SpatOrbs
                    if (tRotateVirtOnly) StartSort = NoOcc + 1
                    if (tRotateOccOnly) EndSort = NoOcc
                else if (x == 2) then
                    if (tSeparateOccVirt) then
                        NoOcc = nOccAlpha
                    else
                        NoOcc = 0
                    end if
                    StartSort = SpatOrbs + 1
                    EndSort = NoOrbs
                    if (tRotateVirtOnly) StartSort = SpatOrbs + NoOcc + 1
                    if (tRotateOccOnly) EndSort = NoOcc + SpatOrbs
                end if

                call sort (Evalues(startSort:endSort), natOrbMat(startSort:endSort, startSort:endSort), &
                           SymOrbs_rotTemp(startSort:endSort))

            end do
               
        else

            ! If we are not truncating, the orbitals get put back into their
            ! original order, so the symmetry information is still correct,
            ! no need for the SymOrbs_rot array. Instead, just take the labels
            ! of SymLabelList3_rot with them.
            do x = 1, NoSpinCyc

                if (x == 1) then
                    if (tSeparateOccVirt) then
                        if (tStoreSpinOrbs) then
                            NoOcc = nOccBeta
                        else
                            NoOcc = NEl/2
                        end if
                    else
                        NoOcc = 0
                    end if
                    StartSort = 1
                    EndSort = SpatOrbs
                    if (tRotateOccOnly) EndSort = NoOcc
                    if (tRotateVirtOnly) StartSort = NoOcc + 1

                else if (x == 2) then
                    if (tSeparateOccVirt) then
                        NoOcc = nOccAlpha
                    else
                        NoOcc =  0
                    end if
                    StartSort = SpatOrbs + 1
                    EndSort = NoOrbs
                    if (tRotateOccOnly) EndSort = SpatOrbs + NoOcc
                    if (tRotateVirtOnly) StartSort = SpatOrbs + NoOcc + 1
                end if

                call sort (EValues(startSort:endSort), NatOrbMat(startSort:endSort, startSort:endSort), &
                           SymLabelList3_rot(startSort:endSort))
            end do 
            
        end if

        call halt_timer(OrderCoeff_Time)

        write(6,*) 'Eigen-values: '
        do i = 1, NoOrbs
            write(6,*) Evalues(i)
        end do

    end subroutine OrderCoeffT1

    subroutine FillCoeffT1

        use RotateOrbsData, only: CoeffT1, SymLabelList3_rot, SymOrbs_rot, SymOrbs_rotTag, &
                                    TruncEval, NoRotOrbs, EvaluesTrunc, EvaluesTruncTag
        use LoggingData, only: tTruncRODump, tTruncDumpbyVal

        integer :: l, k, i, j, NoRotAlphBet, io1, io2, ierr
        character(len=*), parameter :: t_r = 'FillCoeffT1'
        character(len=5) :: Label
        character(len=20) :: LabelFull
        real(dp) :: OccEnergies(1:NoRotOrbs)
  
        FillCoeff_Time%timer_name = 'FillCoeff'
        call set_timer(FillCoeff_Time, 30)

        write(6,*) 'NatOrbMat'
        do i = 1, nBasis
            do j = 1, nBasis
                write(6,'(F10.6)', advance='no') NatOrbMat(j,i)
            end do
            write(6,*) ''
        end do

        if (tTruncRODump) then

            if (tTruncDumpbyVal) then
                NoFrozenVirt = 0
                if (tStoreSpinOrbs) then
                    do i = SpatOrbs, 1, -1
                        if (Evalues(i) > TruncEval) exit
                        if (Evalues(i+SpatOrbs) > TruncEval) exit
                        NoFrozenVirt = NoFrozenVirt + 2
                    end do
                    if (NoFrozenVirt >= (NoOrbs-NEl)) call stop_all(t_r, 'Freezing all virtual orbitals.')
                else
                    do i = SpatOrbs, 1, -1
                        if (Evalues(i) > TruncEval) exit
                        NoFrozenVirt = NoFrozenVirt + 1
                    end do
                    if (NoFrozenVirt >= (SpatOrbs-(NEl/2))) call stop_all(t_r, 'Freezing all virtual orbitals.')
                end if
                NoRotOrbs = NoOrbs - NoFrozenVirt
            end if

            allocate(SymOrbs_rot(NoOrbs), stat=ierr)
            call LogMemAlloc('SymOrbs_rot', NoOrbs, 4, t_r, SymOrbs_rotTag, ierr)
            SymOrbs_rot(:) = 0

            allocate(EvaluesTrunc(NoOrbs-NoFrozenVirt), stat=ierr)
            call LogMemAlloc('EvaluesTrunc', NoOrbs-NoFrozenVirt, 4, t_r, EvaluesTruncTag, ierr)
            EvaluesTrunc(:) = 0.0_dp

            if (tStoreSpinOrbs) then
                NoRotAlphBet = SpatOrbs - (NoFrozenVirt/2)
            else 
                NoRotAlphBet = NoOrbs - NoFrozenVirt
            end if


            if (tStoreSpinOrbs) then                                            
                ! As we reorder these so that they are truncated, we also need
                ! to pair up symmetries.

                call CalcOccEnergies(OccEnergies)

                ! First nOccBeta, then nOccAlpha.
                do i = 1,(2*nOccBeta),2
                    k = 1
                    do while(OccEnergies(k) == 0.0_dp)
                        k = k+2
                    end do
                    do j = 1,(2*nOccBeta),2
                        if ((OccEnergies(j) < OccEnergies(k)) .and. (OccEnergies(j) /= 0.0_dp)) k = j
                    end do
                    l = ceiling(real(k,dp)/2.0_dp)
                    CoeffT1(:,i) = NatOrbMat(:,l)
                    EvaluesTrunc(i) = Evalues(l)
                    SymOrbs_rot(i) = SymOrbs_rotTemp(l)
                    OccEnergies(k) = 0.0_dp
                end do
                do i = 2, (2*nOccAlpha), 2
                    k = 2
                    do while(OccEnergies(k) == 0.0_dp)
                        k = k+2
                    end do
                    do j =2,(2*nOccAlpha),2
                        if ((OccEnergies(j) < OccEnergies(k)).and.(OccEnergies(j) /= 0.0_dp)) k = j
                    end do
                    l= (k/2)+SpatOrbs
                    CoeffT1(:,i) = NatOrbMat(:,l)
                    EvaluesTrunc(i) = Evalues(l)
                    SymOrbs_rot(i) = SymOrbs_rotTemp(l)
                    OccEnergies(k) = 0.0_dp
                end do
                
                ! Need to fill coeffT1 so that it goes alpha beta alpha beta.
                k = (2*nOccBeta)+1
                do i = nOccBeta+1, NoRotAlphBet
                    CoeffT1(:,k) = NatOrbMat(:,i)
                    EvaluesTrunc(k) = Evalues(i)
                    SymOrbs_rot(k) = SymOrbs_rotTemp(i)
                    k = k+2
                end do
                k = (2*nOccAlpha)+2
                do i = SpatOrbs+nOccAlpha+1, SpatOrbs+NoRotAlphBet
                    CoeffT1(:,k) = NatOrbMat(:,i)
                    SymOrbs_rot(k) = SymOrbs_rotTemp(i)
                    EvaluesTrunc(k) = Evalues(i)
                    k = k+2
                end do
 
            else

                ! Order occupied in terms of energy again - this makes sure
                ! freezing etc doesn't get screwed up.
                call CalcOccEnergies(OccEnergies)

                ! OccEnergies has the orbital energies as they are ordered
                ! currently - need to put NatOrbMat into CoeffT1 so that this
                ! goes from lowest energy to highest. 

                do i = 1, NEl/2
                    k = 1
                    do while(OccEnergies(k) == 0.0_dp)
                        k = k+1
                    end do
                    do j = 1, NEl/2
                        if ((OccEnergies(j) < OccEnergies(k)) .and. (OccEnergies(j) /= 0.0_dp)) k = j
                    end do
                    CoeffT1(:,i) = NatOrbMat(:,k)
                    EvaluesTrunc(i) = Evalues(k)
                    SymOrbs_rot(i) = SymOrbs_rotTemp(k)
                    OccEnergies(k) = 0.0_dp
                end do

                do i = (NEl/2)+1, NoRotAlphBet
                    CoeffT1(:,i) = NatOrbMat(:,i)
                    EvaluesTrunc(i) = Evalues(i)
                    SymOrbs_rot(i) = SymOrbs_rotTemp(i)
                end do
            end if

        else

            do i = 1, NoOrbs
                CoeffT1(:,i) = NatOrbMat(:,i)
            end do

        end if

        if (tTruncRODump) then

            Label = ''
            LabelFull = ''
            write(Label,'(I5)') NoFrozenVirt
            LabelFull = 'EVALUES-TRUNC-'//adjustl(Label)

            io1 = get_free_unit()
            open(io1, file=LabelFull, status='unknown')
            if (tStoreSpinOrbs) then
                write(io1,*) NoOrbs-NoFrozenVirt
                do i = 1, NoOrbs-NoFrozenVirt,2
                    write(io1, '(I5,ES20.10,I5,A5,I5,ES20.10,I5)') i, EvaluesTrunc(i), SymOrbs_rot(i), &
                        '  *  ', i+1, EvaluesTrunc(i+1), SymOrbs_rot(i+1)
                end do
            else
                write(io1,*) NoOrbs-NoFrozenVirt
                do i = 1, NoOrbs-NoFrozenVirt
                    write(io1,'(ES20.10,I5)') EvaluesTrunc(i), SymOrbs_rot(i)
                end do
            end if
            close(io1) 
        else
            io2 = get_free_unit()
            open(io2, file='EVALUES', status='unknown')
            write(io2,*) NoOrbs
            if (tStoreSpinOrbs) then
                k = 0
                do i = 1, NoOrbs, 2
                    k = k+1
                    if (tTruncRODump) then
                        write(io2,'(2I5,ES20.10,I5,A5,I5,ES20.10,I5)') (NoOrbs-i+1), i, Evalues(k), SymOrbs_rot(i), &
                            '  *  ', i+1, Evalues(k+SpatOrbs), SymOrbs_rot(i+1)
                    else
                        write(io2,'(2I5,ES20.10,I5,A5,I5,ES20.10,I5)') (NoOrbs-i+1), i, Evalues(k), &
                                int(G1(SymLabelList3_rot(k))%Sym%S,4), '  *  ', &
                             i+1, Evalues(k+SpatOrbs), int(G1(SymLabelList3_rot(k+SpatOrbs))%Sym%S,4)
                    end if
                end do
            else
                do i = 1, SpatOrbs
                    write(io2,'(3I5,ES20.10)') i, NoOrbs-i+1, (NoOrbs-i+1)*2, Evalues(i)
                end do
            end if
            close(io2)
        end if

        call HistNatOrbEvalues()

        call halt_timer(FillCoeff_Time)

    end subroutine FillCoeffT1

    subroutine HistNatOrbEvalues()

        use LoggingData, only: tTruncRODump
        use RotateOrbsData, only: CoeffT1, EvaluesTrunc

        integer :: i, k, a, b, NoOcc, io1, io2
        real(dp) :: EvalueEnergies(1:NoOrbs), OrbEnergies(1:NoOrbs)
        real(dp) :: SumEvalues
        
        io1 = get_free_unit()
        NoOcc = NEl/2 ! Is this correct in all cases?!

        open(io1, file='EVALUES-PLOTRAT', status='unknown')
        if (tStoreSpinOrbs) then
            k = 0
            do i = 1, SpatOrbs
                k = k+2
                write(io1,'(F20.10,ES20.10)') real(k-1,dp)/real(NoOrbs,dp), Evalues(i)
                write(io1,'(F20.10,ES20.10)') real(k,dp)/real(NoOrbs,dp), Evalues(SpatOrbs+i)
            end do
        else if (tRotateOccOnly) then
            k = 0
            do i = 1, NoOcc
                k = k+1
                write(io1,'(F20.10,ES20.10)') real(k,dp)/real(NoOcc,dp), Evalues(i)
            end do
        else if (tRotateVirtOnly) then
            k = NoOcc
            do i = NoOcc+1, NoOrbs
                k = k+1
                write(io1,'(F20.10,ES20.10)') real(k-NoOcc,dp)/real(NoOrbs-NoOcc,dp), Evalues(i)
            end do
        else
            k = 0
            do i = 1, SpatOrbs
                k = k+1
                write(io1,'(F20.10,ES20.10)') real(k,dp)/real(NoOrbs,dp), Evalues(i)
            end do
        end if
        close(io1)

        ! Want to write out the eigenvectors in order of the energy of the new
        ! orbitals - so that we can see the occupations of the type of orbital.
        ! For now, keep this separate to the transformation of ARR - even
        ! though it is equivalent.

        OrbEnergies(:) = 0.0_dp
        EvalueEnergies(:) = 0.0_dp
        SumEvalues = 0.0_dp
        do i = 1, NoOrbs
            if (tStoreSpinOrbs) then
                SumEvalues = SumEvalues + Evalues(i)
            else
                SumEvalues = SumEvalues + (2*Evalues(i))
            end if

            if (tTruncRODump) then
                EvalueEnergies(i) = EvaluesTrunc(i)
            else
                EvalueEnergies(i) = Evalues(i)
            end if

            ! We are only interested in the diagonal elements.
            do a = 1, NoOrbs
                b= SymLabelList2_rot(a)
                if (tStoreSpinOrbs) then
                    OrbEnergies(i) = OrbEnergies(i) + (CoeffT1(a,i)*ARR(b,2)*CoeffT1(a,i))
                else
                    OrbEnergies(i) = OrbEnergies(i) + (CoeffT1(a,i)*ARR(2*b,2)*CoeffT1(a,i))
                end if
            end do
        end do

        ! ARR(:,1) - ordered by energy, ARR(:,2) - ordered by spin-orbital
        ! index.

        call sort (orbEnergies(1:noOrbs), EvalueEnergies(1:noOrbs))

        io2 = get_free_unit()
        open(io2, file='EVALUES-ENERGY', status='unknown')
        do i = 1, NoOrbs
            write(io2,*) OrbEnergies(NoOrbs-i+1), EvalueEnergies(NoOrbs-i+1)
        end do
        write(io2,*) 'The sum of the occupation numbers (eigenvalues) = ', SumEvalues
        write(io2,*) 'The number of electrons = ', NEl
        call neci_flush(io2)
        close(io2)
        call neci_flush(6)

        call PrintOccTable()

    end subroutine HistNatOrbEvalues

    subroutine CalcOccEnergies(OccEnergies)

        use RotateOrbsData, only: CoeffT1, NoRotOrbs

        real(dp) :: OccEnergies(1:NoRotOrbs)
        integer :: i, a, b, NoOcc, x, Prev, k

        OccEnergies(:) = 0.0_dp
        if (tStoreSpinOrbs) then
            do x = 1, 2
                if (x == 1) then
                    NoOcc = nOccBeta
                    Prev = 0
                    k = 1
                else
                    NoOcc = nOccAlpha
                    Prev= SpatOrbs
                    k = 2
                end if
                do i = 1+Prev, NoOcc+Prev
                    ! We are only interested in the diagonal elements.            
                    do a = 1, NoOrbs
                        b= SymLabelList2_rot(a)
                        OccEnergies(k) = OccEnergies(k) + (NatOrbMat(a,i)*ARR(b,2)*NatOrbMat(a,i))
                    end do
                    k = k+2
                end do
            end do
        else
            NoOcc = NEl/2
            do i = 1, NoOcc
                ! We are only interested in the diagonal elements.            
                do a = 1, NoOrbs
                    b = SymLabelList2_rot(a)
                    OccEnergies(i) = OccEnergies(i) + (NatOrbMat(a,i)*ARR(2*b,2)*NatOrbMat(a,i))
                end do
            end do
        end if

    end subroutine CalcOccEnergies

    subroutine PrintOccTable()

        use LoggingData, only: tTruncRODump
        use RotateOrbsData, only: CoeffT1, EvaluesTrunc
        use SystemData, only: tUseHFOrbs

        integer :: x, i, a, b, io2

        io2 = get_free_unit()
        open(io2, file='OccupationTable', status='unknown')
        x = 1
        do while (x <= NoOrbs)
            write(io2,'(A16,A5)', advance='no') 'HF Orb En    ','Sym'
            if (.not.tUseHFOrbs) then
                do i = x, x+9
                    if (i > NoOrbs) then
                        write(io2,*) ''
                        exit
                    end if
                    if (tTruncRODump) then
                        write(io2,'(ES16.6)',advance='no') EvaluesTrunc(i)
                    else
                        write(io2,'(ES16.6)',advance='no') Evalues(i)
                    end if
                end do
            end if
            write(io2,*) ''

            do a = 1, NoOrbs
                b = SymLabelListInv_rot(a)
                if (tStoreSpinOrbs) then
                    write(io2,'(F16.10,I5)',advance='no') ARR(a,1), int(G1(a)%sym%S,4)
                else
                    write(io2,'(F16.10,I5)',advance='no') ARR(2*a,1), int(G1(2*a)%sym%S,4)
                end if
                do i = x, x+9
                    if (i > NoOrbs) then
                        write(io2,*) ''
                        exit
                    end if
                    write(io2,'(F16.10)',advance='no') CoeffT1(b,i)
                end do
                write(io2,*) ''
            end do
            write(io2,*) ''
            x = x+10
        end do
        call neci_flush(io2)
        close(io2)
        call neci_flush(6)
        
    end subroutine PrintOccTable

    subroutine PrintOrbOccs(OrbOccs)

        ! This routine takes whatever orbital basis we're using and is called
        ! at the end of a spawn to find the contribution of each orbital to the
        ! final wavefunction. This is done by histogramming the determinant
        ! populations, and then running over these adding the  coefficients of
        ! each determinant to the orbitals occupied. This is essentially
        ! < Psi | a_p+ a_p | Psi > - the diagonal terms of the one electron
        ! reduced density matrix.

        real(dp) :: Norm,OrbOccs(nBasis),AllOrbOccs(nBasis)
        integer :: i, error, iunit
        logical :: tWarning

        AllOrbOccs = 0.0_dp

        call MPIReduce(OrbOccs,MPI_SUM,AllOrbOccs)

        ! Want to normalise the orbital contributions for convenience.
        tWarning = .false.
        if (iProcIndex == 0) then
            Norm = 0.0_dp
            do i = 1,nBasis
                Norm = Norm + AllOrbOccs(i)
                if ((AllOrbOccs(i) < 0) .or. (Norm < 0)) then
                    write(6,*) 'WARNING: Integer overflow when calculating the orbital occupations.'
                    tWarning = .true.
                end if
            end do
            if (Norm /= 0.0_dp) then
                do i = 1,nBasis
                    AllOrbOccs(i) = AllOrbOccs(i)/Norm
                end do
            end if

            iunit = get_free_unit()
            open(iunit, file = 'ORBOCCUPATIONS', status='UNKNOWN')
            write(iunit, '(A15,A30)') '# Orbital no.','Normalised occupation'
            if (tWarning) write(iunit,*) 'WARNING: integer OVERFLOW OCCURRED WHEN CALCULATING THESE OCCUPATIONS'
            do i = 1,nBasis
                write(iunit,'(I15,F30.10)') i, AllOrbOccs(i)
            end do
            close(iunit)
        end if

    end subroutine PrintOrbOccs

    subroutine DeallocateNatOrbs()

        use LoggingData, only: tTruncRODump

        character(len=*), parameter :: t_r = 'DeallocateNatOrbs'

        if (tTruncRODump) then
            deallocate(SymOrbs_rotTemp)
            call LogMemDeAlloc(t_r, SymOrbs_rotTempTag)
        end if
        deallocate(NatOrbMat)
        call LogMemDeAlloc(t_r, NatOrbMatTag)
        deallocate(Evalues)
        call LogMemDeAlloc(t_r, EvaluesTag)

    end subroutine DeallocateNatOrbs

end module NatOrbsMod
