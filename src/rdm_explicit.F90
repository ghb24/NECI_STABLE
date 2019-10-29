#include "macros.h"

module rdm_explicit

    ! Routines used for calculating the RDM for a pair of FCIQMC wave functions
    ! explicitly, i.e., the RDM is calculated exactly for given wave functions.
    ! Note that this does not mean that the calculated RDMs are exact, since
    ! the wave functions themselves are not exact.

    ! This is slow and should only be used for testing.

    use bit_rep_data, only: NIfTot
    use constants
    use util_mod, only: near_zero
    use SystemData, only : tReltvy, t_3_body_excits, tGUGA, nel
    use guga_bitRepOps, only: encode_matrix_element
    use guga_rdm, only: gen_exc_djs_guga, send_proc_ex_djs, t_test_diagonal
    use bit_reps, only: extract_bit_rep, decode_bit_det

    implicit none

contains

    subroutine Fill_ExplicitRDM_this_Iter(TotWalkers)

        use FciMCData, only: CurrentDets
        use global_utilities, only: set_timer, halt_timer
        use Parallel_neci, only: iProcIndex, MPIAllReduceDataType
        use ParallelHelper, only: MPI_MAXLOC, MPI_2integer
        use rdm_data, only: nElRDM_Time

        integer(int64), intent(in) :: TotWalkers
        integer(n_int) :: iLutnI(0:NIfTot)
        integer(int64) :: MaxTotWalkers, TotWalkIn(2), TotWalkOut(2)
        integer :: i, error
        real(dp) :: TempTotParts, NormalisationTemp, Sum_Coeffs
        logical :: blank_det
        integer :: SignI(lenof_sign), SignI2(lenof_sign)

        ! Run through the current determinants.
        ! Find the max number of determinants on a processor - all need to
        ! run through this number so that the communication can be done at
        ! all stages.

        TotWalkIn(1) = TotWalkers
        TotWalkIn(2) = iProcIndex

        call MPIAllReduceDatatype(TotWalkIn, 1, MPI_MAXLOC, MPI_2integer, TotWalkOut)

        MaxTotWalkers = TotWalkOut(1)

        call set_timer(nElRDM_Time, 30)

        do i = 1, int(MaxTotWalkers,sizeof_int)

            ! But if the actual number of determinants on this processor is
            ! less than the number  we're running through, feed in 0
            ! determinants and 0 sign.
            if (i .gt. TotWalkers) then
                iLutnI(:) = 0
                blank_det = .true.
            else
                iLutnI(:) = CurrentDets(:,i)
                blank_det = .false.
            end if

            call Add_ExplicitRDM_Contrib(iLutnI, blank_det)

        end do

        call halt_timer(nElRDM_Time)

    end subroutine Fill_ExplicitRDM_this_Iter

    subroutine Fill_Hist_ExplicitRDM_this_Iter()

        use bit_reps, only: encode_sign
        use DetCalcData, only: Det, FCIDets
        use hist_data, only: AllHistogram, Histogram
        use global_utilities, only: set_timer, halt_timer
        use Parallel_neci, only: iProcIndex, MPISumAll
        use rdm_data, only: nElRDM_Time, ExcNorm

        integer(n_int) :: iLutnI(0:NIfTot)
        integer :: i, error
        real(dp) :: TempTotParts, NormalisationTemp, Sum_Coeffs, AllNode_norm
        logical :: blank_det
        real(dp), dimension(lenof_sign) :: TempSign

        call set_timer(nElRDM_Time,30)

        call MPISumAll(Histogram, AllHistogram)

        ExcNorm = 0.0_dp
        if (iProcIndex .eq. 0) then
            do i = 1, Det
                ExcNorm = ExcNorm + AllHistogram(1,i)**2
            end do
            ExcNorm = sqrt(ExcNorm)
        end if

        call MPISumAll(ExcNorm, allNode_norm)
        ExcNorm = allNode_norm

        do i = 1, Det

            ! But if the actual number of determinants on this processor is
            ! less than the number we're running through, feed in 0
            ! determinants and 0 sign.
            if (near_zero(Histogram(1, i))) then
                iLutnI(:) = 0
                blank_det = .true.
            else
                iLutnI(:) = FCIDets(:,i)
                blank_det = .false.
            end if

            TempSign(1) = real(i,dp)
            call encode_sign(iLutnI, TempSign)

            call Add_Hist_ExplicitRDM_Contrib(iLutnI, blank_det)

        end do

        call halt_timer(nElRDM_Time)

    end subroutine Fill_Hist_ExplicitRDM_this_Iter

    subroutine Add_ExplicitRDM_Contrib(iLutnI, blank_det)

        ! This is the general routine for taking a particular determinant in
        ! the spawned list, D_i and adding it's contribution to the reduced
        ! density matrix.

        use LoggingData, only: RDMExcitLevel
        use Parallel_neci, only: nProcessors
        use rdm_data, only: Sing_ExcList, Doub_ExcList, Sing_InitExcSlots, Doub_InitExcSlots
        use rdm_data, only: Sing_ExcDjs, Doub_ExcDjs
        use bit_rep_data, only: nifguga
        use guga_bitRepOps, only: convert_ilut_toGUGA

        integer(n_int), intent(in) :: iLutnI(0:NIfTot)
        logical, intent(in) :: blank_det
        integer :: i, ni(nel), FlagsDi
        integer(n_int) :: ilutG(0:nifguga)
        real(dp) :: SignDi(lenof_sign)

        ! Set up excitation arrays.
        ! These are blocked according to the processor the excitation would be
        ! on if occupied. In each block, the first entry is the sign of
        ! determinant D_i and the second the bit string of the determinant
        ! (these need to be sent along with the excitations). Each processor
        ! will have a different Di.

        Sing_ExcDjs(:,:) = 0
        Sing_ExcList(:) = 0
        Sing_ExcList(:) = Sing_InitExcSlots(:)

        if (tGUGA) then
            call convert_ilut_toGUGA(ilutNi, ilutG)

            call extract_bit_rep(iLutnI, nI, SignDi, FlagsDi)
            call encode_matrix_element(ilutG, SignDi(1), 2)

            do i = 0, nProcessors - 1
                Sing_ExcDjs(:, Sing_ExcList(i)) = ilutG
                Sing_ExcList(i) = Sing_ExcList(i) + 1
            end do

        else

            do i = 0, nProcessors-1
                Sing_ExcDjs(:,Sing_ExcList(i)) = iLutnI(:)
                Sing_ExcList(i) = Sing_ExcList(i) + 1
            end do
        end if

        if (RDMExcitLevel /= 1) then
            Doub_ExcDjs(:,:) = 0
            Doub_ExcList(:) = 0
            Doub_ExcList(:) = Doub_InitExcSlots(:)

            if (tGUGA) then
                do i = 0, nProcessors - 1
                    Doub_ExcDjs(:,Doub_ExcList(i)) = ilutG
                    Doub_ExcList(i) = Doub_ExcList(i) + 1
                end do
            else
                do i = 0,nProcessors-1
                    Doub_ExcDjs(:,Doub_ExcList(i)) = iLutnI(:)
                    Doub_ExcList(i) = Doub_ExcList(i) + 1
                end do
            end if
        end if

        ! Out of here we will get a filled ExcDjs array with all the single or
        ! double excitations from Dj, this will be done for each proc.
        if (.not. blank_det) then
            if (tGUGA) then
                call gen_exc_djs_guga(ilutnI)
            else
                call GenExcDjs(iLutnI)
            end if
        end if

        ! We then need to send the excitations to the relevant processors.
        ! This routine then calls SearchOccDets which takes each excitation
        ! and and binary searches the occupied determinants for this. If found,
        ! we re-find the orbitals and parity involved in the excitation, and
        ! add the c_i*c_j contributions to the corresponding matrix element.
        if (tGUGA) then
            call send_proc_ex_djs()
        else
            call SendProcExcDjs()
        end if

    end subroutine Add_ExplicitRDM_Contrib

    subroutine Add_Hist_ExplicitRDM_Contrib(iLutnI,blank_det)

        ! This is the general routine for taking a particular determinant in
        ! the spawned list, D_i and adding it's contribution to the reduced
        ! density matrix.

        use LoggingData, only: RDMExcitLevel
        use Parallel_neci, only: nProcessors
        use rdm_data, only: Sing_ExcList, Doub_ExcList, Sing_InitExcSlots, Doub_InitExcSlots
        use rdm_data, only: Sing_ExcDjs, Doub_ExcDjs

        integer(n_int), intent(in) :: iLutnI(0:NIfTot)
        logical, intent(in) :: blank_det
        integer :: i

        ! Set up excitation arrays.
        ! These are blocked according to the processor the excitation would be
        ! on if occupied. In each block, the first entry is the sign of
        ! determinant D_i and the second the bit string of the determinant
        ! (these need to be sent along with the excitations). Each processor
        ! will have a different Di.

        Sing_ExcDjs(:,:) = 0
        Sing_ExcList(:) = 0
        Sing_ExcList(:) = Sing_InitExcSlots(:)

        do i = 0, nProcessors-1
            Sing_ExcDjs(:,Sing_ExcList(i)) = iLutnI(:)
            Sing_ExcList(i) = Sing_ExcList(i)+1
        end do
        if (RDMExcitLevel /= 1) then
            Doub_ExcDjs(:,:) = 0
            Doub_ExcList(:) = 0
            Doub_ExcList(:) = Doub_InitExcSlots(:)

            do i = 0, nProcessors-1
                Doub_ExcDjs(:,Doub_ExcList(i)) = iLutnI(:)
                Doub_ExcList(i) = Doub_ExcList(i)+1
            end do
        end if

        ! Out of here we will get a filled ExcDjs array with all the single or
        ! double excitations  from Dj, this will be done for each proc.
        if (.not. blank_det) call Gen_Hist_ExcDjs(iLutnI)

        ! We then need to send the excitations to the relevant processors.
        ! This routine then calls SearchOccDets which takes each excitation
        ! and and binary searches the occupied determinants for this. If found,
        ! we re-find the orbitals and parity involved in the excitation, and add
        ! the c_i*c_j contributions to the corresponding matrix element.
        call Send_Hist_ProcExcDjs()

    end subroutine Add_Hist_ExplicitRDM_Contrib

    subroutine GenExcDjs(iLutnI)

        ! This uses GenExcitations3 in symexcit3.F90 to generate all the
        ! possible either single or double excitations from D_i, finds the
        ! processor they would be on if occupied, and puts them in the
        ! SingExcDjs array according to that processor.

        use DetBitOps, only: EncodeBitDet
        use load_balance_calcnodes, only: DetermineDetNode
        use LoggingData, only: RDMExcitLevel
        use Parallel_neci, only: nProcessors
        use rdm_data, only: Sing_ExcList, Doub_ExcList, OneEl_Gap, TwoEl_Gap
        use rdm_data, only: Sing_ExcDjs, Doub_ExcDjs
        use rdm_data, only: one_rdms, two_rdm_spawn
        use rdm_filling, only: fill_diag_1rdm, fill_spawn_rdm_diag
        use SymExcit3, only: GenExcitations3
        use SymExcit4, only: GenExcitations4, ExcitGenSessionType
        use SystemData, only: nel, tReltvy

        integer(n_int), intent(in) :: iLutnI(0:NIfTot)

        integer(n_int) :: iLutnJ(0:NIfTot)
        real(dp) :: SignDi(lenof_sign), full_sign(1)
        integer :: ExcitMat3(2,2), nI(NEl), nJ(NEl), Proc, FlagsDi
        integer :: a, b, exflag
        logical :: tAllExcitFound, tParity

        type(ExcitGenSessionType) :: session

        ! Unfortunately uses the decoded determinant - might want to look at this.
        call extract_bit_rep(iLutnI, nI, SignDi, FlagsDi)

        if (RDMExcitLevel == 1) then
            call fill_diag_1rdm(one_rdms, nI, SignDi)
        else
            full_sign = SignDi(1)*SignDi(lenof_sign)
            call fill_spawn_rdm_diag(two_rdm_spawn, nI, full_sign)
        end if

        if (.not. t_test_diagonal) then
        ! Zeros in ExcitMat3 starts off at the first single excitation.
        ExcitMat3(:,:) = 0

        ! This becomes true when all the excitations have been found.
        tAllExcitFound = .false.

        do while (.not. tAllExcitFound)
            exflag = 1

            ! Passed out of here is the singly excited determinant, nJ.
            ! Information such as the orbitals involved in the excitation and
            ! the parity is also found in this step, we are not currently
            ! storing this, and it is re-calculated later on (after the
            ! determinants are passed to the relevant processor) - but the
            ! speed of sending this information vs recalculating it will be
            ! tested. RDMExcitLevel is passed through, if this is 1, only
            ! singles are generated, if it is 2 only doubles are found.
            if (tReltvy) then
                call GenExcitations4(session, nI, nJ, exFlag, ExcitMat3(:,:), tParity, tAllExcitFound, .true.)
            else
                call GenExcitations3(nI, iLutnI, nJ, exflag, ExcitMat3(:,:), tParity, tAllExcitFound, .true.)
            endif

            if (tAllExcitFound) exit

            iLutnJ(:) = 0
            call EncodeBitDet(nJ, iLutnJ)

            Proc = DetermineDetNode(nel, nJ, 0)
            ! This will return a value between 0 -> nProcessors-1
            Sing_ExcDjs(:,Sing_ExcList(Proc)) = iLutnJ(:)
            Sing_ExcList(Proc) = Sing_ExcList(Proc)+1

            ! Want a quick test to see if arrays are getting full.
            if (Sing_ExcList(Proc) .gt. nint(OneEl_Gap*(Proc+1))) then
                write(6,*) 'Proc', Proc
                write(6,*) 'Sing_ExcList', Sing_ExcList
                write(6,*) 'No. spaces for each proc', nint(OneEl_Gap)
                call Stop_All('GenExcDjs', 'Too many excitations for space available.')
            end if
        end do

        if (RDMExcitLevel /= 1) then

            ! Zeros in ExcitMat3 starts off at the first single excitation.
            ExcitMat3(:,:) = 0

            ! This becomes true when all the excitations have been found.
            tAllExcitFound = .false.

            do while (.not. tAllExcitFound)
                exflag = 2

                ! Passed out of here is the doubly excited determinant, nJ.
                ! Information such as the orbitals involved in the excitation
                ! and the parity is  also found in this step, we are not
                ! currently storing this, and it is re-calculated later on
                ! (after the determinants are passed to the relevant processor)
                ! - but the speed of sending this information vs recalculating
                ! it will be tested. RDMExcitLevel is passed through, if this
                ! is 1, only singles are generated, if it is 2 only doubles are
                if (tReltvy) then
                    call GenExcitations4(session, nI, nJ, exFlag, ExcitMat3(:,:), tParity, tAllExcitFound, .true.)
                else
                    ! found.
                    call GenExcitations3(nI, iLutnI, nJ, exflag, ExcitMat3(:,:), tParity, tAllExcitFound, .true.)
                endif

                if (tAllExcitFound) exit

                iLutnJ(:) = 0
                call EncodeBitDet(nJ, iLutnJ)

                Proc = DetermineDetNode(nel, nJ, 0)
                !This will return a value between 0 -> nProcessors-1
                Doub_ExcDjs(:,Doub_ExcList(Proc)) = iLutnJ(:)
                ! All the double excitations from this particular nI are being
                ! stored in Doub_ExcDjs.

                Doub_ExcList(Proc) = Doub_ExcList(Proc)+1

                ! Want a quick test to see if arrays are getting full.
                if (Doub_ExcList(Proc) .gt. nint(TwoEl_Gap*(Proc+1))) then
                    write(6,*) 'Proc', Proc
                    write(6,*) 'Doub_ExcList', Doub_ExcList
                    write(6,*) 'No. spaces for each proc', nint(TwoEl_Gap)
                    call Stop_All('GenExcDjs','Too many excitations for space available.')
                end if
            end do
        end if
        end if

    end subroutine GenExcDjs

    subroutine Gen_Hist_ExcDjs(iLutnI)

        ! This uses GenExcitations3 in symexcit3.F90 to generate all the
        ! possible either single or double excitations from D_i, finds the
        ! processor they would be on if occupied, and puts them in the
        ! SingExcDjs array according to that processor.

        use DetBitOps, only: EncodeBitDet
        use hist_data, only: AllHistogram
        use load_balance_calcnodes, only: DetermineDetNode
        use LoggingData, only: RDMExcitLevel
        use Parallel_neci, only: nProcessors
        use rdm_data, only: Sing_ExcList, Doub_ExcList, ExcNorm, OneEl_Gap, TwoEl_Gap
        use rdm_data, only: Sing_ExcDjs, Doub_ExcDjs
        use rdm_data, only: one_rdms, two_rdm_spawn
        use rdm_filling, only: fill_diag_1rdm, fill_spawn_rdm_diag
        use SymExcit3, only: GenExcitations3
        use SymExcit4, only: GenExcitations4, ExcitGenSessionType
        use SystemData, only: nel

        integer(n_int), intent(in) :: iLutnI(0:NIfTot)
        integer(n_int) :: iLutnJ(0:NIfTot)
        integer, dimension(lenof_sign) :: HistPos
        real(dp), dimension(lenof_sign) :: RealHistPos
        integer :: ExcitMat3(2,2), nI(NEl), nJ(NEl), Proc, FlagsDi
        integer :: a, b, exflag
        logical :: tAllExcitFound, tParity
        real(dp) :: SignDi(lenof_sign), full_sign(1)

        type(ExcitGenSessionType) :: session

        ! Unfortunately uses the decoded determinant - might want to look at this.
        call extract_bit_rep (iLutnI, nI, RealHistPos, FlagsDi)

        HistPos = int(RealHistPos)

        SignDi(1) = AllHistogram(1, HistPos(1))/ExcNorm
        SignDi(lenof_sign) = AllHistogram(1,HistPos(1))/ExcNorm

        if (RDMExcitLevel == 1) then
            call fill_diag_1rdm(one_rdms, nI, SignDi)
        else
            full_sign = SignDi(1)*SignDi(lenof_sign)
            call fill_spawn_rdm_diag(two_rdm_spawn, nI, full_sign)
        end if

        ! Zeros in ExcitMat3 starts off at the first single excitation.
        ExcitMat3(:,:) = 0

        ! This becomes true when all the excitations have been found.
        tAllExcitFound = .false.

        do while (.not. tAllExcitFound)
            exflag = 1

            ! Passed out of here is the singly excited determinant, nJ.
            ! Information such as the orbitals involved in the excitation and
            ! the parity is also found in this step, we are not currently
            ! storing this, and it is re-calculated later on (after the
            ! determinants are passed to the relevant processor) - but the speed
            ! of sending this information vs recalculating it will be tested.
            ! RDMExcitLevel is passed through, if this is 1, only singles are
            ! generated, if it is 2 only doubles are found.
            if (tReltvy) then
                call GenExcitations4(session, nI, nJ, exFlag, ExcitMat3(:,:), tParity, tAllExcitFound, .true.)
            else
                call GenExcitations3(nI, iLutnI, nJ, exflag, ExcitMat3(:,:), tParity, tAllExcitFound, .true.)
            endif

            if (tAllExcitFound) exit

            iLutnJ(:) = 0
            call EncodeBitDet(nJ,iLutnJ)

            Proc = DetermineDetNode(nel,nJ,0)
            ! This will return a value between 0 -> nProcessors-1.
            Sing_ExcDjs(:,Sing_ExcList(Proc)) = iLutnJ(:)
            Sing_ExcList(Proc) = Sing_ExcList(Proc) + 1

            ! Want a quick test to see if arrays are getting full.
            if (Sing_ExcList(Proc) .gt. nint(OneEl_Gap*(Proc+1))) then
                write(6,*) 'Proc', Proc
                write(6,*) 'Sing_ExcList', Sing_ExcList
                write(6,*) 'No. spaces for each proc', nint(OneEl_Gap)
                call Stop_All('GenExcDjs', 'Too many excitations for space available.')
            end if
        end do

        if (RDMExcitLevel /= 1) then

            ExcitMat3(:,:) = 0
            ! Zeros in ExcitMat3 starts off at the first single excitation.
            tAllExcitFound = .false.
            ! This becomes true when all the excitations have been found.

            do while (.not. tAllExcitFound)
                exflag = 2

                ! Passed out of here is the doubly excited determinant, nJ.
                ! Information such as the orbitals involved in the excitation
                ! and the parity is also found in this step, we are not currently
                ! storing this, and it is re-calculated  later on (after the
                ! determinants are passed to the relevant processor) - but the
                ! speed of sending this information vs recalculating it will be
                ! tested. RDMExcitLevel is passed through, if this is 1, only
                ! singles are generated, if it is 2 only doubles are found.
                 if (tReltvy) then
                    call GenExcitations4(session, nI, nJ, exFlag, ExcitMat3(:,:), tParity, tAllExcitFound, .true.)
                else
                    call GenExcitations3(nI, iLutnI, nJ, exflag, ExcitMat3(:,:), tParity, tAllExcitFound, .true.)
                endif

                if (tAllExcitFound) exit

                iLutnJ(:) = 0
                call EncodeBitDet(nJ,iLutnJ)

                Proc = DetermineDetNode(nel, nJ, 0)
                ! This will return a value between 0 -> nProcessors-1.
                Doub_ExcDjs(:,Doub_ExcList(Proc)) = iLutnJ(:)
                ! All the double excitations from this particular nI are being
                ! stored in Doub_ExcDjs.

                Doub_ExcList(Proc) = Doub_ExcList(Proc) + 1

                ! Want a quick test to see if arrays are getting full.
                if (Doub_ExcList(Proc) .gt. nint(TwoEl_Gap*(Proc+1))) then
                    write(6,*) 'Proc', Proc
                    write(6,*) 'Doub_ExcList', Doub_ExcList
                    write(6,*) 'No. spaces for each proc', nint(TwoEl_Gap)
                    call Stop_All('GenExcDjs','Too many excitations for space available.')
                end if
            end do
        end if

    end subroutine Gen_Hist_ExcDjs

    subroutine SendProcExcDjs()

        ! In this routine the excitations are sent to the relevant processors.
        ! Sent with them will be the Di they were excited from and its sign.
        ! Each processor will receive nProcessor number of lists with different
        ! Di determinants. The original Di's will (I think) still be in the
        ! original InitSingExcSlots positions. This follows the
        ! directannihilation algorithm closely.

        use LoggingData, only: RDMExcitLevel
        use Parallel_neci, only: nProcessors, MPIArg, MPIAlltoAll, MPIAlltoAllv
        use rdm_data, only: Sing_ExcList, Doub_ExcList, OneEl_Gap, TwoEl_Gap
        use rdm_data, only: Sing_ExcDjs, Doub_ExcDjs, Sing_ExcDjs2, Doub_ExcDjs2

        integer :: i, j
        integer :: error, MaxSendIndex,MaxIndex
        integer(MPIArg) :: sendcounts(nProcessors),disps(nProcessors)
        integer(MPIArg) :: sing_recvcounts(nProcessors)
        integer(MPIArg) :: sing_recvdisps(nProcessors)
        integer(MPIArg) :: doub_recvcounts(nProcessors),doub_recvdisps(nProcessors)

        do i = 0, nProcessors - 1
            ! Sendcounts is the number of singly excited determinants we want
            ! to send for each processor (but goes from 1, not 0).
            sendcounts(i+1) = int(Sing_ExcList(i)-(nint(OneEl_Gap*i)+1), MPIArg)

            ! disps is the first slot for each processor - 1.
            disps(i+1) = nint(OneEl_Gap*i, MPIArg)
        end do

        MaxSendIndex = Sing_ExcList(nProcessors-1) - 1

        ! We now need to calculate the recvcounts and recvdisps -
        ! this is a job for AlltoAll
        sing_recvcounts(1:nProcessors) = 0
        call MPIAlltoAll(sendcounts, 1, sing_recvcounts, 1, error)

        ! We can now get recvdisps from recvcounts, since we want the data to
        ! be contiguous after the move.
        sing_recvdisps(1) = 0
        do i = 2, nProcessors
            sing_recvdisps(i) = sing_recvdisps(i-1) + sing_recvcounts(i-1)
        end do

        MaxIndex = sing_recvdisps(nProcessors) + sing_recvcounts(nProcessors)
        ! But the actual number of integers we need to send is the calculated
        ! values * NIfTot+1.
        do i = 1, nProcessors
            sendcounts(i) = sendcounts(i)*(int(NIfTot+1,MPIArg))
            disps(i) = disps(i)*(int(NIfTot+1,MPIArg))
            sing_recvcounts(i) = sing_recvcounts(i)*(int(NIfTot+1,MPIArg))
            sing_recvdisps(i) = sing_recvdisps(i)*(int(NIfTot+1,MPIArg))
        end do

#ifdef PARALLEL
        call MPIAlltoAllv(Sing_ExcDjs(:,1:MaxSendIndex), sendcounts, disps,&
                            Sing_ExcDjs2, sing_recvcounts, sing_recvdisps, error)
#else
        Sing_ExcDjs2(0:NIfTot,1:MaxIndex) = Sing_ExcDjs(0:NIfTot,1:MaxSendIndex)
#endif

        call Sing_SearchOccDets(sing_recvcounts, sing_recvdisps)

        if (RDMExcitLevel /= 1) then
            do i = 0, nProcessors-1
                ! Sendcounts is the number of singly excited determinants we
                ! want to send for each processor (but goes from 1, not 0).
                sendcounts(i+1) = int(Doub_ExcList(i)-(nint(TwoEl_Gap*i)+1),MPIArg)

                ! disps is the first slot for each processor - 1.
                disps(i+1) = nint(TwoEl_Gap*i,MPIArg)
            end do

            MaxSendIndex = Doub_ExcList(nProcessors-1) - 1

            ! We now need to calculate the recvcounts and recvdisps -
            ! this is a job for AlltoAll
            doub_recvcounts(1:nProcessors) = 0
            call MPIAlltoAll(sendcounts, 1, doub_recvcounts, 1, error)

            ! We can now get recvdisps from recvcounts, since we want the data
            ! to be contiguous after the move.
            doub_recvdisps(1) = 0
            do i = 2, nProcessors
                doub_recvdisps(i) = doub_recvdisps(i-1) + doub_recvcounts(i-1)
            end do

            MaxIndex = doub_recvdisps(nProcessors)+doub_recvcounts(nProcessors)
            ! But the actual number of integers we need to send is the
            ! calculated values * NIfTot+1.
            do i = 1, nProcessors
                sendcounts(i) = sendcounts(i)*(int(NIfTot+1,MPIArg))
                disps(i) = disps(i)*(int(NIfTot+1,MPIArg))
                doub_recvcounts(i) = doub_recvcounts(i)*(int(NIfTot+1,MPIArg))
                doub_recvdisps(i) = doub_recvdisps(i)*(int(NIfTot+1,MPIArg))
            end do

            ! This is the main send of all the single excitations to the
            ! corresponding processors.
#ifdef PARALLEL
            call MPIAlltoAllv(Doub_ExcDjs(:,1:MaxSendIndex), sendcounts, disps,&
                                    Doub_ExcDjs2, doub_recvcounts, doub_recvdisps, error)
#else
            Doub_ExcDjs2(0:NIfTot,1:MaxIndex) = Doub_ExcDjs(0:NIfTot,1:MaxSendIndex)
#endif
            call Doub_SearchOccDets(doub_recvcounts, doub_recvdisps)

        end if

    end subroutine SendProcExcDjs

    subroutine Send_Hist_ProcExcDjs()

        ! In this routine the excitations are sent to the relevant processors.
        ! Sent with them will be the Di they were excited from and its sign.
        ! Each processor will receive nProcessor number of lists with different
        ! Di determinants. The original Di's will (I think) still be in the
        ! original InitSingExcSlots positions. This follows the
        ! directannihilation algorithm closely.

        use LoggingData, only: RDMExcitLevel
        use Parallel_neci, only: nProcessors, MPIArg, MPIAlltoAll, MPIAlltoAllv
        use rdm_data, only: Sing_ExcList, Doub_ExcList, OneEl_Gap, TwoEl_Gap
        use rdm_data, only: Sing_ExcDjs, Doub_ExcDjs, Sing_ExcDjs2, Doub_ExcDjs2

        integer :: i, j
        integer :: error, MaxSendIndex, MaxIndex
        integer(MPIArg) :: sendcounts(nProcessors),disps(nProcessors)
        integer(MPIArg) :: sing_recvcounts(nProcessors)
        integer(MPIArg) :: sing_recvdisps(nProcessors)
        integer(MPIArg) :: doub_recvcounts(nProcessors),doub_recvdisps(nProcessors)

        do i = 0, nProcessors - 1
            ! sendcounts is the number of singly excited determinants we want to send for
            ! each processor (but goes from 1, not 0).
            sendcounts(i+1) = int(Sing_ExcList(i)-(nint(OneEl_Gap*i)+1), MPIArg)

            ! disps is the first slot for each processor - 1.
            disps(i+1) = nint(OneEl_Gap*i, MPIArg)
        end do

        MaxSendIndex = Sing_ExcList(nProcessors-1) - 1

        ! We now need to calculate the recvcounts and recvdisps -
        ! this is a job for AlltoAll.
        sing_recvcounts(1:nProcessors)=0
        call MPIAlltoAll(sendcounts, 1, sing_recvcounts, 1, error)

        ! We can now get recvdisps from recvcounts, since we want the data to
        ! be contiguous after the move.
        sing_recvdisps(1) = 0
        do i = 2, nProcessors
            sing_recvdisps(i) = sing_recvdisps(i-1)+sing_recvcounts(i-1)
        end do

        MaxIndex = sing_recvdisps(nProcessors)+sing_recvcounts(nProcessors)
        ! But the actual number of integers we need to send is the calculated
        ! values * NIfTot+1.
        do i = 1, nProcessors
            sendcounts(i) = sendcounts(i)*(int(NIfTot+1,MPIArg))
            disps(i) = disps(i)*(int(NIfTot+1,MPIArg))
            sing_recvcounts(i) = sing_recvcounts(i)*(int(NIfTot+1,MPIArg))
            sing_recvdisps(i) = sing_recvdisps(i)*(int(NIfTot+1,MPIArg))
        end do
#ifdef PARALLEL
        call MPIAlltoAllv(Sing_ExcDjs(:,1:MaxSendIndex), sendcounts, disps, Sing_ExcDjs2, sing_recvcounts, sing_recvdisps, error)
#else
        Sing_ExcDjs2(0:NIfTot,1:MaxIndex) = Sing_ExcDjs(0:NIfTot,1:MaxSendIndex)
#endif

        call Sing_Hist_SearchOccDets(sing_recvcounts, sing_recvdisps)


        if (RDMExcitLevel /= 1) then
            do i = 0, nProcessors-1
                ! sendcounts is the number of singly excited determinants we
                ! want to send for each processor (but goes from 1, not 0).
                sendcounts(i+1) = int(Doub_ExcList(i)-(nint(TwoEl_Gap*i)+1),MPIArg)

                ! disps is the first slot for each processor - 1.
                disps(i+1) = nint(TwoEl_Gap*i,MPIArg)
            end do

            MaxSendIndex = Doub_ExcList(nProcessors-1)-1

            ! We now need to calculate the recvcounts and recvdisps -
            ! this is a job for AlltoAll.
            doub_recvcounts(1:nProcessors) = 0
            call MPIAlltoAll(sendcounts, 1, doub_recvcounts, 1, error)

            ! We can now get recvdisps from recvcounts, since we want the data
            ! to be contiguous after the move.
            doub_recvdisps(1) = 0
            do i = 2, nProcessors
                doub_recvdisps(i) = doub_recvdisps(i-1) + doub_recvcounts(i-1)
            end do

            MaxIndex = doub_recvdisps(nProcessors) + doub_recvcounts(nProcessors)
            ! But the actual number of integers we need to send is the calculated
            ! values * NIfTot+1.
            do i = 1, nProcessors
                sendcounts(i) = sendcounts(i)*(int(NIfTot+1,MPIArg))
                disps(i) = disps(i)*(int(NIfTot+1,MPIArg))
                doub_recvcounts(i) = doub_recvcounts(i)*(int(NIfTot+1,MPIArg))
                doub_recvdisps(i) = doub_recvdisps(i)*(int(NIfTot+1,MPIArg))
            end do

            ! This is the main send of all the single excitations to the
            ! corresponding processors.
#ifdef PARALLEL
            call MPIAlltoAllv(Doub_ExcDjs(:,1:MaxSendIndex),sendcounts,disps,&
                                    Doub_ExcDjs2,doub_recvcounts,doub_recvdisps,error)
#else
            Doub_ExcDjs2(0:NIfTot,1:MaxIndex)=Doub_ExcDjs(0:NIfTot,1:MaxSendIndex)
#endif
            call Doub_Hist_SearchOccDets(doub_recvcounts,doub_recvdisps)

        end if

    end subroutine Send_Hist_ProcExcDjs

    subroutine Sing_SearchOccDets(recvcounts,recvdisps)

        ! We now have arrays SingExcDjs2 which contain all the single
        ! excitations from each processor. These number sent from processor i
        ! is recvcounts(i), and the first 2 have information about the
        ! determinant Di from which the Dj's are single excitations (and it's sign).

        use FciMCData, only: CurrentDets, TotWalkers
        use LoggingData, only: RDMExcitLevel
        use Parallel_neci, only: nProcessors, MPIArg
        use rdm_data, only: Sing_ExcDjs, Sing_ExcDjs2
        use rdm_data, only: one_rdms, two_rdm_spawn
        use rdm_filling, only: fill_sings_1rdm, fill_spawn_rdm_singles
        use searching, only: BinSearchParts_rdm
        use SystemData, only: nel
#ifdef __DEBUG
        character(*), parameter :: this_routine = "Sing_SearchOccDets"
#endif

        integer(MPIArg), intent(in) :: recvcounts(nProcessors),recvdisps(nProcessors)

        integer(n_int) :: iLutnJ(0:NIfTot)
        real(dp) :: SignDi(lenof_sign), SignDj(lenof_sign), full_sign(1)
        integer :: i, j, NoDets, StartDets, PartInd
        integer :: nI(NEl), nJ(NEl), Ex(2,2), Ex_symm(2,2), FlagsDi, FlagsDj
        logical :: tDetFound, tParity

        ! Take each Dj, and binary search CurrentDets to see if it is occupied.
        do i = 1, nProcessors
            ! Doing determinants from each processor separately because each
            ! has a different D_i it was excited from.
            NoDets = recvcounts(i)/(NIfTot+1)
            StartDets = (recvdisps(i)/(NIfTot+1))+1

            if (NoDets .gt. 1) then
                call extract_bit_rep(Sing_ExcDjs2(:,StartDets), nI, SignDi, FlagsDi)

                do j = StartDets+1, (NoDets+StartDets-1)

                    ! D_i is in the first spot - start from the second.
                    iLutnJ(:) = Sing_ExcDjs2(:,j)
                    ! This binary searches CurrentDets between 1 and
                    ! TotWalkers for determinant iLutnJ. If found, tDetFound
                    ! will be true, and PartInd the index in CurrentDets where
                    ! the determinant is.
                    call BinSearchParts_rdm(iLutnJ,1,int(TotWalkers,sizeof_int),PartInd,tDetFound)
                    if (tDetFound) then
                        ! Determinant occupied; add c_i*c_j to the relevant
                        ! element of nElRDM. Need to first find the orbitals
                        ! involved in the excitation from D_i -> D_j and the parity.
                        Ex(:,:) = 0
                        ! Ex(1,1) goes in as the max number of excitations -
                        ! we know this is an excitation of level RDMExcitLevel.
                        Ex(1,1) = 1
                        tParity = .false.

                        call extract_bit_rep(CurrentDets(:,PartInd), nJ, SignDj, FlagsDj)

                        ! Ex(1,:) comes out as the orbital(s) excited from,
                        ! Ex(2,:) comes out as the orbital(s) excited to.
                        call GetExcitation(nI, nJ, NEl, Ex, tParity)

                        if (Ex(1,1) .le. 0) call Stop_All('Sing_SearchOccDets', 'nJ is not the correct excitation of nI.')

                        if (RDMExcitLevel == 1) then
                            call fill_sings_1rdm(one_rdms, Ex, tParity, SignDi, SignDj, .true.)
                        else
                            ASSERT(.not. t_3_body_excits)
                            if (tParity) then
                                full_sign = -SignDi(1)*SignDj(lenof_sign)
                            else
                                full_sign = SignDi(1)*SignDj(lenof_sign)
                            end if

                            call fill_spawn_rdm_singles(two_rdm_spawn, nI, Ex, full_sign)
                            ! Add in symmetric contribution.
                            Ex_symm(1,:) = Ex(2,:)
                            Ex_symm(2,:) = Ex(1,:)
                            call fill_spawn_rdm_singles(two_rdm_spawn, nI, Ex_symm, full_sign)
                        end if

                        ! No normalisation factor just yet - possibly need to revise.
                    end if

                end do
            end if
        end do

    end subroutine Sing_SearchOccDets

    subroutine Doub_SearchOccDets(recvcounts,recvdisps)

        ! We now have arrays SingExcDjs2 which contain all the single excitations
        ! from each processor. These number sent from processor i is
        ! recvcounts(i), and the first 2 have information  about the determinant
        ! Di from which the Dj's are single excitations (and it's sign).

        use FciMCData, only: CurrentDets, TotWalkers
        use Parallel_neci, only: nProcessors, MPIArg
        use rdm_data, only: Doub_ExcDjs, Doub_ExcDjs2, two_rdm_spawn
        use rdm_data_utils, only: add_to_rdm_spawn_t
        use searching, only: BinSearchParts_rdm
        use SystemData, only: nel
#ifdef __DEBUG
        character(*), parameter :: this_routine = "Doub_SearchOccDets"
#endif
        integer(MPIArg), intent(in) :: recvcounts(nProcessors),recvdisps(nProcessors)
        integer(n_int) :: iLutnJ(0:NIfTot)
        real(dp) :: SignDi(lenof_sign), SignDj(lenof_sign), full_sign(1)
        integer :: i, j, NoDets, StartDets, PartInd
        integer :: nI(NEl), nJ(NEl), Ex(2,2), FlagsDi, FlagsDj
        logical :: tDetFound, tParity

        ! Take each Dj, and binary search CurrentDets to see if it is occupied.
        do i = 1, nProcessors
            ! Doing determinants from each processor separately because each
            ! has a different D_i it was excited from.
            NoDets = recvcounts(i)/(NIfTot+1)
            StartDets = (recvdisps(i)/(NIfTot+1))+1

            if (NoDets.gt.1) then
                call extract_bit_rep (Doub_ExcDjs2(:,StartDets), nI, SignDi, FlagsDi)

                do j = StartDets+1, (NoDets+StartDets-1)
                    ! D_i is in the first spot - start from the second.
                    iLutnJ(:) = Doub_ExcDjs2(:,j)

                    ! This binary searches CurrentDets between 1 and TotWalkers
                    ! for determinant iLutnJ. If found, tDetFound will be true,
                    ! and PartInd the index in CurrentDets where the determinant is.
                    call BinSearchParts_rdm(iLutnJ, 1, int(TotWalkers,sizeof_int), PartInd, tDetFound)

                    if (tDetFound) then
                        ! Determinant occupied; add c_i*c_j to the relevant
                        ! element of nElRDM. Need to first find the orbitals
                        ! involved in the excitation from D_i -> D_j and
                        ! the parity.
                        Ex(:,:) = 0
                        ! Ex(1,1) goes in as the max number of excitations -
                        ! we know this is an excitation of level RDMExcitLevel.
                        Ex(1,1) = 2
                        tParity = .false.

                        call extract_bit_rep (CurrentDets(:,PartInd), nJ, SignDj, FlagsDj)

                        ! Ex(1,:) comes out as the orbital(s) excited from,
                        ! Ex(2,:) comes out as the orbital(s) excited to.
                        call GetExcitation(nI, nJ, NEl, Ex, tParity)

                        if (Ex(1,1) .le. 0) call Stop_All('SearchOccDets',&
                                            'nJ is not the correct excitation of nI.')

                        ASSERT(.not. t_3_body_excits)
                        if (tParity) then
                            full_sign = -SignDi(1)*SignDj(lenof_sign)
                        else
                            full_sign = SignDi(1)*SignDj(lenof_sign)
                        end if

                        call add_to_rdm_spawn_t(two_rdm_spawn, Ex(2,1), Ex(2,2), Ex(1,1), Ex(1,2), full_sign, .false.)
                        ! Add in symmetric contribution.
                        call add_to_rdm_spawn_t(two_rdm_spawn, Ex(1,1), Ex(1,2), Ex(2,1), Ex(2,2), full_sign, .false.)
                    end if
                end do
            end if
        end do

    end subroutine Doub_SearchOccDets

    subroutine Sing_Hist_SearchOccDets(recvcounts,recvdisps)

        ! We now have arrays SingExcDjs2 which contain all the single
        ! excitations from each processor. These number sent from processor i
        ! is recvcounts(i), and the first 2 have information about the
        ! determinant Di from which the Dj's are single excitations (and it's sign).

        use DetBitOps, only: FindBitExcitLevel
        use FciMCData, only: ilutHF_True, TotWalkers
        use hist, only: find_hist_coeff_explicit
        use hist_data, only: AllHistogram
        use LoggingData, only: RDMExcitLevel
        use Parallel_neci, only: nProcessors, MPIArg
        use rdm_data, only: Sing_ExcDjs, Sing_ExcDjs2, ExcNorm
        use rdm_data, only: one_rdms, two_rdm_spawn
        use rdm_filling, only: fill_sings_1rdm, fill_spawn_rdm_singles
        use searching, only: BinSearchParts_rdm
        use SystemData, only: nel

        integer(MPIArg), intent(in) :: recvcounts(nProcessors),recvdisps(nProcessors)
        integer(n_int) :: iLutnJ(0:NIfTot)
        integer, dimension(lenof_sign) :: HistPos
        real(dp), dimension(lenof_sign) :: RealHistPos

        integer :: i, j, NoDets, StartDets, PartInd, ExcitLevel
        integer :: nI(NEl), nJ(NEl), Ex(2,2), Ex_symm(2,2), FlagsDi, FlagsDj
        logical :: tDetFound, tParity
        real(dp) :: SignDi(lenof_sign), SignDj(lenof_sign), full_sign(1)

#ifdef __DEBUG
        character(*), parameter :: this_routine = "Sing_Hist_SearchOccDets"
#endif
        ! Take each Dj, and binary search CurrentDets to see if it is occupied.
        do i = 1, nProcessors

            ! Doing determinants from each processor separately because each
            ! has a different D_i it was excited from.
            NoDets = recvcounts(i)/(NIfTot+1)
            StartDets = (recvdisps(i)/(NIfTot+1)) + 1

            if (NoDets .gt. 1) then
                call extract_bit_rep (Sing_ExcDjs2(:,StartDets), nI, RealHistPos, FlagsDi)

                HistPos=int(RealHistPos)

                SignDi = AllHistogram(1,HistPos(1))/ExcNorm

                do j=StartDets+1,(NoDets+StartDets-1)

                    ! D_i is in the first spot - start from the second.
                    iLutnJ(:) = Sing_ExcDjs2(:,j)
                    ! This binary searches CurrentDets between 1 and TotWalkers
                    ! for determinant iLutnJ. If found, tDetFound will be true,
                    ! and PartInd the index in CurrentDets where the determinant is.
                    call BinSearchParts_rdm(iLutnJ,1,int(TotWalkers,sizeof_int),PartInd,tDetFound)

                    ExcitLevel = FindBitExcitLevel (iLutHF_true, iLutnJ, NEl)
                    call find_hist_coeff_explicit (iLutnJ, ExcitLevel, PartInd, tDetFound)

                    if (tDetFound) then
                        ! Determinant occupied; add c_i*c_j to the relevant
                        ! element of nElRDM. Need to first find the orbitals
                        ! involved in the excitation from D_i -> D_j and the
                        ! parity.
                        Ex(:,:) = 0
                        ! Ex(1,1) goes in as the max number of excitations - we
                        ! know this is an excitation of level RDMExcitLevel.
                        Ex(1,1) = 1
                        tParity = .false.

                        call decode_bit_det(nJ,iLutnJ)

                        SignDj = AllHistogram(1,PartInd)/ExcNorm

                        ! Ex(1,:) comes out as the orbital(s) excited from,
                        ! Ex(2,:) comes out as the orbital(s) excited to.
                        call GetExcitation(nI, nJ, NEl, Ex, tParity)

                        if (Ex(1,1).le.0) call Stop_All('Sing_SearchOccDets',&
                                            'nJ is not the correct excitation of nI.')

                        if (RDMExcitLevel == 1) then
                            call fill_sings_1rdm(one_rdms, Ex, tParity, SignDi, SignDj, .true.)
                        else
                            ASSERT(.not. t_3_body_excits)
                            if (tParity) then
                                full_sign = -SignDi(1)*SignDj(lenof_sign)
                            else
                                full_sign = SignDi(1)*SignDj(lenof_sign)
                            end if

                            call fill_spawn_rdm_singles(two_rdm_spawn, nI, Ex, full_sign)
                            ! Add in symmetric contribution.
                            Ex_symm(1,:) = Ex(2,:)
                            Ex_symm(2,:) = Ex(1,:)
                            call fill_spawn_rdm_singles(two_rdm_spawn, nI, Ex_symm, full_sign)
                        end if

                        ! No normalisation factor just yet - possibly need to revise.
                    end if

                end do
            end if
        end do

    end subroutine Sing_Hist_SearchOccDets

    subroutine Doub_Hist_SearchOccDets(recvcounts,recvdisps)

        ! We now have arrays SingExcDjs2 which contain all the single excitations
        ! from each processor. These number sent from processor i is recvcounts(i),
        ! and the first 2 have information about the determinant Di from which
        ! the Dj's are single excitations (and it's sign).

        use DetBitOps, only: FindBitExcitLevel
        use FciMCData, only: ilutHF_True, TotWalkers
        use hist, only: find_hist_coeff_explicit
        use hist_data, only: AllHistogram
        use Parallel_neci, only: nProcessors, MPIArg
        use rdm_data, only: Doub_ExcDjs, Doub_ExcDjs2, ExcNorm, two_rdm_spawn
        use rdm_data_utils, only: add_to_rdm_spawn_t
        use searching, only: BinSearchParts_rdm
        use SystemData, only: nel

        integer(MPIArg), intent(in) :: recvcounts(nProcessors),recvdisps(nProcessors)
        integer(n_int) :: iLutnJ(0:NIfTot)
        integer, dimension(lenof_sign) :: HistPos
        real(dp), dimension(lenof_sign) :: RealHistPos
        integer :: i, j, NoDets, StartDets,PartInd, ExcitLevel
        integer :: nI(NEl), nJ(NEl), Ex(2,2), FlagsDi, FlagsDj
        logical :: tDetFound, tParity
        real(dp) :: SignDi(lenof_sign), SignDj(lenof_sign), full_sign(1)

#ifdef __DEBUG
        character(*), parameter :: this_routine = "Doub_Hist_SearchOccDets"
#endif
        ! Take each Dj, and binary search the CurrentDets to see if it is occupied.
        do i = 1, nProcessors

            ! Doing determinants from each processor separately because each
            ! has a different D_i it was excited from.
            NoDets=recvcounts(i)/(NIfTot+1)
            StartDets=(recvdisps(i)/(NIfTot+1))+1

            if (NoDets .gt. 1) then
                call extract_bit_rep (Doub_ExcDjs2(:,StartDets), nI, RealHistPos, FlagsDi)

                HistPos = int(RealHistPos)

                SignDi = AllHistogram(1,HistPos(1))/ExcNorm

                do j = StartDets + 1, (NoDets+StartDets-1)

                    ! D_i is in the first spot - start from the second.
                    iLutnJ(:) = Doub_ExcDjs2(:,j)

                    ! This binary searches CurrentDets between 1 and TotWalkers
                    ! for determinant iLutnJ. If found, tDetFound will be true,
                    ! and PartInd the index in CurrentDets where the  determinant is.
                    call BinSearchParts_rdm(iLutnJ, 1, int(TotWalkers,sizeof_int), PartInd, tDetFound)

                    ExcitLevel = FindBitExcitLevel (iLutHF_True, iLutnJ, NEl)
                    call find_hist_coeff_explicit (iLutnJ, ExcitLevel, PartInd, tDetFound)

                    if (tDetFound) then
                        ! Determinant occupied; add c_i*c_j to the relevant
                        ! element of nElRDM. Need to first find the orbitals
                        ! involved in the excitation from D_i -> D_j and
                        ! the parity.
                        Ex(:,:) = 0
                        ! Ex(1,1) goes in as the max number of excitations - we
                        ! know this is an excitation of level RDMExcitLevel.
                        Ex(1,1) = 2
                        tParity = .false.

                        call decode_bit_det(nJ,iLutnJ)
                        SignDj = AllHistogram(1,PartInd)/ExcNorm

                        ! Ex(1,:) comes out as the orbital(s) excited from,
                        ! Ex(2,:) comes out as the orbital(s) excited to.
                        call GetExcitation(nI, nJ, NEl, Ex, tParity)

                        if (Ex(1,1) .le. 0) call Stop_All('SearchOccDets', 'nJ is not the correct excitation of nI.')

                        ASSERT(.not. t_3_body_excits)
                        if (tParity) then
                            full_sign = -SignDi(1)*SignDj(lenof_sign)
                        else
                            full_sign = SignDi(1)*SignDj(lenof_sign)
                        end if

                        call add_to_rdm_spawn_t(two_rdm_spawn, Ex(2,1), Ex(2,2), Ex(1,1), Ex(1,2), full_sign, .false.)
                        ! Add in symmetric contribution.
                        call add_to_rdm_spawn_t(two_rdm_spawn, Ex(1,1), Ex(1,2), Ex(2,1), Ex(2,2), full_sign, .false.)
                    end if
                end do
            end if
        end do

    end subroutine Doub_Hist_SearchOccDets

end module rdm_explicit
