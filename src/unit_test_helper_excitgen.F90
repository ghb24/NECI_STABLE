#include "macros.h"

module unit_test_helper_excitgen
    use constants
    use bit_reps, only: IlutBits, init_bit_rep
    use read_fci, only: readfciint, initfromfcid, fcidump_name
    use shared_memory_mpi, only: shared_allocate_mpi, shared_deallocate_mpi
    use IntegralsData, only: UMat, umat_win
    use Integrals_neci, only: IntInit, get_umat_el_normal
    use procedure_pointers, only: get_umat_el, generate_excitation
    use SystemData, only: nel, nBasis, UMatEps, tStoreSpinOrbs, tReadFreeFormat, &
                          tReadInt, t_pcpp_excitgen
    use sort_mod
    use System, only: SysInit, SetSysDefaults, SysCleanup
    use Parallel_neci, only: MPIInit, MPIEnd
    use UMatCache, only: GetUMatSize, tTransGTID, setupUMat2d_dense
    use OneEInts, only: Tmat2D
    use bit_rep_data, only: NIfTot, nifd, extract_sign
    use bit_reps, only: encode_sign, decode_bit_det
    use DetBitOps, only: EncodeBitDet, DetBitEq
    use SymExcit3, only: countExcitations3, GenExcitations3
    use FciMCData, only: pSingles, pDoubles, pParallel, ilutRef, projEDet, &
                         fcimc_excit_gen_store
    use SymExcitDataMod, only: excit_gen_store_type, scratchSize
    use GenRandSymExcitNUMod, only: init_excit_gen_store, construct_class_counts
    use Calc, only: CalcInit, CalcCleanup, SetCalcDefaults
    use dSFMT_interface, only: dSFMT_init, genrand_real2_dSFMT
    use Determinants, only: DetInit, DetPreFreezeInit, get_helement, DefDet, tDefineDet
    use util_mod, only: get_free_unit
    use orb_idx_mod, only: SpinProj_t
    implicit none

    abstract interface
        function calc_pgen_t(nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2) result(pgen)
            use constants
            use SymExcitDataMod, only: scratchSize
            use bit_rep_data, only: NIfTot
            use SystemData, only: nel
            implicit none
            integer, intent(in) :: nI(nel)
            integer(n_int), intent(in) :: ilutI(0:NIfTot)
            integer, intent(in) :: ex(2, 2), ic
            integer, intent(in) :: ClassCount2(ScratchSize), ClassCountUnocc2(ScratchSize)

            real(dp) :: pgen

        end function calc_pgen_t
    end interface

    procedure(calc_pgen_t), pointer :: calc_pgen

    abstract interface
        subroutine to_unit_writer_t(iunit)
            integer, intent(in) :: iunit
        end subroutine
    end interface

    type, abstract :: Writer_t
        procedure(to_unit_writer_t), pointer, nopass :: write
        ! I would like it to be:
        ! character(:), allocatable :: filepath
        ! but for gfortran <= 4.8.5 it has to be
        character(512) :: filepath
    end type

    type, extends(Writer_t) :: FciDumpWriter_t
    end type

    type, extends(Writer_t) :: InputWriter_t
    end type

contains

    subroutine test_excitation_generator(sampleSize, pTot, pNull, numEx, nFound, t_calc_pgen, start_nI)
        ! Test an excitation generator by creating a large number of excitations and
        ! compare the generated excits with a precomputed list of all excitations
        ! We thus make sure that
        !   a) all possible excitations are generated with some weight
        !   b) no invalid excitations are obtained
        implicit none
        integer, intent(in) :: sampleSize
        real(dp), intent(out) :: pTot, pNull
        integer, intent(out) :: numEx, nFound
        logical, intent(in) :: t_calc_pgen
        integer, intent(in), optional :: start_nI(nEl)
        integer :: nI(nel), nJ(nel)
        integer :: i, ex(2, maxExcit), exflag
        integer(n_int) :: ilut(0:NIfTot), ilutJ(0:NIfTot)
        real(dp) :: pgen
        logical :: tPar, tAllExFound, tFound
        integer :: j, nSingles, nDoubles
        integer(n_int), allocatable :: allEx(:, :)
        real(dp) :: pgenArr(lenof_sign)
        real(dp) :: matel, matelN, pgenCalc
        logical :: exDoneDouble(0:nBasis, 0:nBasis, 0:nBasis, 0:nBasis)
        logical :: exDoneSingle(0:nBasis, 0:nBasis)
        integer :: ic, part, nullExcits
        integer :: ClassCountOcc(scratchSize), ClassCountUnocc(scratchSize)
        integer(int64) :: start, finish, rate
        character(*), parameter :: t_r = "test_excitation_generator"
        HElement_t(dp) :: HEl
        exDoneDouble = .false.
        exDoneSingle = .false.
        call system_clock(count_rate=rate)

        ! some starting det - do NOT use the reference for the pcpp test, that would
        ! defeat the purpose
        if (present(start_nI)) then
            nI = start_nI
        else
            do i = 1, nel
                if (2 * i < nBasis) then
                    nI(i) = 2 * i - mod(i, 2)
                else
                    nI(i) = i
                end if
            end do
            call sort(nI)
        end if

        call EncodeBitDet(nI, ilut)

        exflag = 3
        ex = 0
        ! create a list of all singles and doubles for reference
        call CountExcitations3(nI, exflag, nSingles, nDoubles)
        allocate(allEx(0:(NIfTot + 1), nSingles + nDoubles), source=0_n_int)
        numEx = 0
        tAllExFound = .false.
        do
            call GenExcitations3(nI, ilut, nJ, exflag, ex, tPar, tAllExFound, .false.)
            if (tAllExFound) exit
            call encodeBitDet(nJ, ilutJ)
            numEx = numEx + 1
            allEx(0:nifd, numEx) = ilutJ(0:nifd)
        end do
        call sort(nI)

        write(iout, *) "In total", numEx, "excits, (", nSingles, nDoubles, ")"
        write(iout, *) "Exciting from", nI

        call EncodeBitDet(nI, ilut)

        ! set the biases for excitation generation
        pParallel = 0.5_dp
        pSingles = 0.1_dp
        pDoubles = 0.9_dp
        pNull = 0.0_dp
        nullExcits = 0
        call system_clock(start)
        do i = 1, sampleSize
            fcimc_excit_gen_store%tFilled = .false.
            call generate_excitation(nI, ilut, nJ, ilutJ, exFlag, ic, ex, tPar, pgen, HEl, fcimc_excit_gen_store, part)
            ! lookup the excitation
            tFound = .false.
            do j = 1, numEx
                if (DetBitEQ(ilutJ, allEx(:, j))) then
                    pgenArr = pgen
                    call encode_sign(allEx(:, j), pgenArr)
                    ! increase its counter by 1
                    allEx(NIfTot + 1, j) = allEx(NIfTot + 1, j) + 1
                    tFound = .true.
                    exit
                end if
            end do
            ! if it was not found, and is not marked as invalid, we have a problem: this is not
            ! an excitaion
            if (.not. tFound .and. .not. nJ(1) == 0) then
                call decode_bit_det(nJ, ilutJ)
                write(iout, *) "Created excitation", nJ
                call stop_all(t_r, "Error: Invalid excitation")
            end if
            ! check if the generated excitation is invalid, if it is, mark this specific constellation
            ! so we do not double-count when calculating pNull
            if (nJ(1) == 0) then
                nullExcits = nullExcits + 1
                if (ic == 2) then
                    if (.not. exDoneDouble(ex(1, 1), ex(1, 2), ex(2, 1), ex(2, 2))) then
                        exDoneDouble(ex(1, 1), ex(1, 2), ex(2, 1), ex(2, 2)) = .true.
                        pNull = pNull + pgen
                    end if
                else if (ic == 1) then
                    if (.not. exDoneSingle(ex(1, 1), ex(1, 2))) then
                        exDoneSingle(ex(1, 1), ex(1, 2)) = .true.
                        pNull = pNull + pgen
                    end if
                end if
            end if
        end do
        call system_clock(finish)

        ! check that all excits have been generated and all pGens are right
        ! probability normalization
        pTot = pNull
        matelN = 0.0
        do i = 1, numEx
            call decode_bit_det(nJ, allEx(:, i))
            matelN = matelN + abs(get_helement(nI, nJ))
        end do
        nFound = 0
        ! class counts might be required for comparing the pgen
        call construct_class_counts(nI, classCountOcc, classCountUnocc)
        do i = 1, numEx
            call extract_sign(allEx(:, i), pgenArr)
            call decode_bit_det(nJ, allEx(:, i))
            matel = get_helement(nI, nJ)
            if (pgenArr(1) > eps) then
                nFound = nFound + 1
                write(iout, *) i, pgenArr(1), real(allEx(NIfTot + 1, i)) / real(sampleSize), &
                    abs(matel) / (pgenArr(1) * matelN)
                ! compare the stored pgen to the directly computed one
                if (t_calc_pgen) then
                    if (i > nSingles) then
                        ic = 2
                    else
                        ic = 1
                    end if
                    ex(1, 1) = 2
                    call getBitExcitation(ilut, allEx(:, i), ex, tPar)
                    pgenCalc = calc_pgen(nI, ilut, ex, ic, ClassCountOcc, ClassCountUnocc)
                    if (abs(pgenArr(1) - pgenCalc) > eps) then
                        write(iout, *) "Stored: ", pgenArr(1), "calculated:", pgenCalc
                        write(iout, *) "For excit", nJ
                        call stop_all(t_r, "Incorrect pgen")
                    end if
                end if
            else
                ! excitations with zero matrix element are not required to be found
                if (abs(matel) < eps) then
                    nFound = nFound + 1
                else if (i < nSingles) then
                    write(iout, *) "Unfound single excitation", nJ
                else
                    write(iout, *) "Unfound double excitation", nJ, matel
                end if
            end if
            pTot = pTot + pgenArr(1)
        end do
        write(iout, *) "Total prob. ", pTot
        write(iout, *) "pNull ", pNull
        write(iout, *) "Null ratio", nullExcits / real(sampleSize)
        write(iout, *) "In total", numEx, "excitations"
        write(iout, *) "With", nSingles, "single excitation"
        write(iout, *) "Found", nFound, "excitations"
        write(iout, *) 'Elapsed Time in seconds:', dble(finish - start) / dble(rate)
        write(iout, *) 'Elapsed Time in micro seconds per excitation:', dble(finish - start) * 1e6_dp / dble(sampleSize* rate)

    end subroutine test_excitation_generator

    !------------------------------------------------------------------------------------------!

    subroutine init_excitgen_test(ref_det, fcidump_writer)
        ! mimick the initialization of an FCIQMC calculation to the point where we can generate
        ! excitations with a weighted excitgen
        ! This requires setup of the basis, the symmetries and the integrals
        integer, intent(in) :: ref_det(:)
        type(FciDumpWriter_t), intent(in) :: fcidump_writer
        integer :: nBasisMax(5, 3), lms
        integer(int64) :: umatsize
        real(dp) :: ecore
        character(*), parameter :: this_routine = 'init_excitgen_test'
        integer, parameter :: seed = 25

        umatsize = 0
        nel = size(ref_det)

        IlutBits%len_orb = 0
        IlutBits%ind_pop = 1
        IlutBits%len_pop = 1
        IlutBits%len_tot = 2

        nifd = 0
        NIfTot = 2

        fcidump_name = "FCIDUMP"
        UMatEps = 1.0e-8
        tStoreSpinOrbs = .false.
        tTransGTID = .false.
        tReadFreeFormat = .true.

        call dSFMT_init(seed)

        call SetCalcDefaults()
        call SetSysDefaults()
        tReadInt = .true.

        call write_file(fcidump_writer)

        get_umat_el => get_umat_el_normal

        call initfromfcid(nel, nbasismax, nBasis, lms, .false.)

        call GetUMatSize(nBasis, umatsize)

        allocate(TMat2d(nBasis, nBasis))

        call shared_allocate_mpi(umat_win, umat, [umatsize])

        UMat = h_cast(0._dp)
        call readfciint(UMat, umat_win, nBasis, ecore, .false.)

        ! init the umat2d storage
        call setupUMat2d_dense(nBasis)

        call SysInit()
        ! required: set up the spin info

        call DetInit()
        ! call SpinOrbSymSetup()

        tDefineDet = .true.
        DefDet = ref_det
        call DetPreFreezeInit()

        call CalcInit()

        call set_ref()

        t_pcpp_excitgen = .true.
        call init_excit_gen_store(fcimc_excit_gen_store)

        call init_bit_rep()
    end subroutine init_excitgen_test

    !------------------------------------------------------------------------------------------!

    subroutine finalize_excitgen_test()
        deallocate(TMat2D)
        call shared_deallocate_mpi(umat_win, UMat)
        call CalcCleanup()
        call SysCleanup()
    end subroutine finalize_excitgen_test

    !------------------------------------------------------------------------------------------!

    ! generate an FCIDUMP file with random numbers with a given sparsity and write to iunit
    subroutine generate_random_integrals(iunit, n_el, n_spat_orb, sparse, sparseT, total_ms)
        integer, intent(in) :: iunit, n_el, n_spat_orb
        real(dp), intent(in) :: sparse, sparseT
        type(SpinProj_t), intent(in) :: total_ms
        integer :: i, j, k, l
        real(dp) :: r, matel
        ! we get random matrix elements from the cauchy-schwartz inequalities, so
        ! only <ij|ij> are random -> random 2d matrix
        real(dp) :: umatRand(n_spat_orb, n_spat_orb)

        umatRand = 0.0_dp
        do i = 1, n_spat_orb
            do j = 1, n_spat_orb
                r = genrand_real2_dSFMT()
                if (r < sparse) &
                    umatRand(i, j) = r * r
            end do
        end do

        ! write the canonical FCIDUMP header
        write(iunit, *) "&FCI NORB=", n_spat_orb, ",NELEC=", n_el, "MS2=", total_ms%val, ","
        write(iunit, "(A)", advance="no") "ORBSYM="
        do i = 1, n_spat_orb
            write(iunit, "(A)", advance="no") "1,"
        end do
        write(iunit, *)
        write(iunit, *) "ISYM=1,"
        write(iunit, *) "&END"
        ! generate random 4-index integrals with a given sparsity
        do i = 1, n_spat_orb
            do j = 1, i
                do k = i, n_spat_orb
                    do l = 1, k
                        matel = sqrt(umatRand(i, j) * umatRand(k, l))
                        if (matel > eps) write(iunit, *) matel, i, j, k, l
                    end do
                end do
            end do
        end do
        ! generate random 2-index integrals with a given sparsity
        do i = 1, n_spat_orb
            do j = 1, i
                r = genrand_real2_dSFMT()
                if (r < sparseT) then
                    write(iunit, *) r, i, j, 0, 0
                end if
            end do
        end do
    end subroutine generate_random_integrals

    !------------------------------------------------------------------------------------------!
    subroutine generate_uniform_integrals

        use SystemData, only: nSpatOrbs, nel, lms

        integer :: i, j, k, l, iunit

        iunit = get_free_unit()

        open (iunit, file="FCIDUMP")
        write(iunit, *) "&FCI NORB=", nSpatOrbs, ",NELEC=", nel, "MS2=", lms, ","
        write(iunit, "(A)", advance="no") "ORBSYM="
        do i = 1, nSpatOrbs
            write(iunit, "(A)", advance="no") "1,"
        end do
        write(iunit, *)
        write(iunit, *) "ISYM=1,"
        write(iunit, *) "&END"

        do i = 1, nSpatOrbs
            do j = 1, nSpatOrbs
                do l = 1, nSpatOrbs
                    do k = 1, nSpatOrbs
                        write(iunit, *) h_cast(1.0_dp), i, j, k, l
                    end do
                end do
            end do
        end do
        do i = 1, nSpatOrbs
            do j = i, nSpatOrbs
                write(iunit, *) h_cast(1.0_dp), i, j, 0, 0
            end do
        end do

        write(iunit, *) h_cast(0.0_dp), 0, 0, 0, 0

        close (iunit)

    end subroutine generate_uniform_integrals

    !------------------------------------------------------------------------------------------!

    ! set the reference to the determinant with the first nel orbitals occupied
    subroutine set_ref()
        integer :: i

        projEDet = reshape([(i + 2, i = 1, nel)], [nel, 1])
        if (allocated(ilutRef)) deallocate(ilutRef)
        allocate(ilutRef(0:NifTot, 1))
        call encodeBitDet(projEDet(:, 1), ilutRef(:, 1))
    end subroutine

    subroutine free_ref()
        deallocate(ilutRef)
        deallocate(projEDet)
    end subroutine free_ref

    subroutine delete_file(path)
        character(*), intent(in) :: path
        integer :: file_id

        file_id = get_free_unit()
        open (file_id, file=path, status='old')
        close (file_id, status='delete')
    end subroutine

    subroutine write_file(writer)
        class(Writer_t), intent(in) :: writer
        integer :: file_id

        file_id = get_free_unit()
        open (file_id, file=writer%filepath)
        call writer%write(file_id)
        close (file_id)
    end subroutine
end module unit_test_helper_excitgen
