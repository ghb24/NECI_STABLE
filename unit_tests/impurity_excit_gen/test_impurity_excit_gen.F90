#include "macros.h"

module test_impurity_excit_gen_mod
    use constants, only: dp
    use SystemData, only: t_complex_ints, tUHF
    use unit_test_helper_excitgen, only: FciDumpWriter_t

    better_implicit_none
    type, extends(FciDumpWriter_t) :: ImpurityFciDumpWriter_t
        integer :: nBasis, nEl, ms
    contains
        procedure :: write_to_unit
    end type


contains

    subroutine write_to_unit(this, iunit)
        class(ImpurityFciDumpWriter_t), intent(in) :: this
        integer, intent(in) :: iunit
        integer :: i
        real(dp) :: e
        unused_var(this)
        ! write the canonical FCIDUMP header
        tUHF = .true.
        write(iunit, *) "&FCI NORB=", this%nBasis, ",NELEC=", this%nel, "MS2=", this%ms, ","
        write(iunit, "(A)", advance="no") "ORBSYM="
        do i = 1, this%nBasis
            write(iunit, "(A)", advance="no") "1,"
        end do
        write(iunit, *)
        write(iunit, *) "ISYM=1,UHF=T"
        write(iunit, *) "&END"
        write(iunit, *) 1, 1, 1, 2, 2
        do i = 1, this%nBasis/2-1
            write(iunit, *) 0.1, 1, 2*i+1, 0, 0
            write(iunit, *) 0.1, 2, 2*i+2, 0, 0
            write(iunit, *) 0.1, 2*i+1, 1, 0, 0
            write(iunit, *) 0.1, 2*i+2, 2, 0, 0
        end do
        do i = 1, this%nBasis
            if(i < 3) then
                e = -0.5
            else if(i < this%nBasis/2+2) then
                e = -0.3
            else
                e = 0.3
            endif
            write(iunit, *) e, i, i, 0, 0
        end do
    end subroutine

end module

program test_impurity_excit_gen
    use constants, only: dp
    use impurity_models, only: gen_excit_impurity_model, setupImpurityExcitgen
    use OneEInts, only: tCPMDSymTMat, tOneElecDiag
    use procedure_pointers, only: generate_excitation
    use Parallel_neci, only: MPIInit, MPIEnd
    use fruit, only: init_fruit, fruit_summary, fruit_finalize, &
        get_failed_count, run_test_case, assert_true, assert_equals
    use util_mod, only: stop_all
    use unit_test_helper_excitgen, only: init_excitgen_test, test_excitation_generator
    use SystemData, only: nBasis, t_complex_ints, tUHF

    use test_impurity_excit_gen_mod, only: ImpurityFciDumpWriter_t

    implicit none
    integer, parameter :: nBasis_ = 6
    integer, parameter :: nel_ = 3
    integer, parameter :: ms_ = 1


    integer :: failed_count

    call MPIInit(.false.)
    call init_fruit()
    ! this is where the tests go
    call impurity_excit_gen_test_driver
    call fruit_summary
    call fruit_finalize

    call get_failed_count(failed_count)

    if (failed_count /= 0) call stop_all('test_impurity_excit_gen', 'failed_tests')

    call MPIEnd(.false.)
contains

    subroutine impurity_excit_gen_test_driver
        implicit none
        integer, parameter :: nSamples = 1000
        real(dp) :: pTot, pNull
        integer :: nFound, numEx, i

        nBasis = nBasis_
        tCPMDSymTMat = .false.
        tOneElecDiag = .false.
        t_complex_ints = .false.

        associate(ref_det => [(i, i = 1, nel_)])
            call init_excitgen_test(ref_det, ImpurityFciDumpWriter_t(nBasis=nBasis_, nEl=nEl_, ms=ms_, filepath='FCIDUMP'))
        end associate
        call init_impurity_tests
        call test_excitation_generator(nSamples, pTot, pNull, numEx, nFound, .false.)
        call assert_equals(numEx,nFound)
        call assert_true(abs(1.0-pTot) < 0.001)
    end subroutine impurity_excit_gen_test_driver


    !------------------------------------------------------------------------------------------!

    subroutine init_impurity_tests
      implicit none

      generate_excitation => gen_excit_impurity_model
      call setupImpurityExcitgen()

    end subroutine init_impurity_tests

end program
