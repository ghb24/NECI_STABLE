#include "macros.h"

program test_impurity_excit_gen
    use impurity_models
    use fruit
    use unit_test_helper_excitgen

    implicit none

    integer, parameter :: nBasis_ = 6
    integer, parameter :: nel_ = 3
    integer, parameter :: ms_ = 1

    call MPIInit(.false.)
    call init_fruit()
    ! this is where the tests go
    call impurity_excit_gen_test_driver
    call fruit_summary
    call fruit_finalize
    call MPIEnd(.false.)

contains

    subroutine impurity_excit_gen_test_driver
        implicit none
        integer, parameter :: nSamples = 1000
        real(dp) :: pTot, pNull
        integer :: nFound, numEx, i
        call init_excitgen_test([(i, i = 1, nel_)], FciDumpWriter_t(impurity_fcidump, 'FCIDUMP'))
        call init_impurity_tests
        call test_excitation_generator(nSamples, pTot, pNull, numEx, nFound, .false.)
        call assert_equals(numEx,nFound)
        call assert_true(abs(1.0-pTot) < 0.001)
    end subroutine impurity_excit_gen_test_driver

    !------------------------------------------------------------------------------------------!

    subroutine impurity_fcidump(iunit)
        use OneEInts, only: tCPMDSymTMat, tOneElecDiag
        use SystemData, only: nBasis, t_complex_ints
        integer, intent(in) :: iunit
        integer :: i
        real(dp) :: e
        nBasis = nBasis_
        tCPMDSymTMat = .false.
        tOneElecDiag = .false.
        t_complex_ints = .false.
        ! write the canonical FCIDUMP header
        write(iunit, *) "&FCI NORB=", nBasis, ",NELEC=", nel_, "MS2=", ms_, ","
        write(iunit, "(A)", advance="no") "ORBSYM="
        do i = 1, nBasis
            write(iunit, "(A)", advance="no") "1,"
        end do
        write(iunit, *)
        write(iunit, *) "ISYM=1,UHF=T"
        write(iunit, *) "&END"
        write(iunit, *) 1, 1, 1, 2, 2
        do i = 1, nBasis/2-1
            write(iunit, *) 0.1, 1, 2*i+1, 0, 0
            write(iunit, *) 0.1, 2, 2*i+2, 0, 0
            write(iunit, *) 0.1, 2*i+1, 1, 0, 0
            write(iunit, *) 0.1, 2*i+2, 2, 0, 0
        end do
        do i = 1, nBasis
            if(i < 3) then
                e = -0.5
            else if(i < nBasis/2+2) then
                e = -0.3
            else
                e = 0.3
            endif
            write(iunit, *) e, i, i, 0, 0
        end do
    end subroutine impurity_fcidump

    !------------------------------------------------------------------------------------------!

    subroutine init_impurity_tests
      implicit none

      generate_excitation => gen_excit_impurity_model
      call setupImpurityExcitgen()

    end subroutine init_impurity_tests

    ! subroutine generate_test_ilut(ilut,nI)
    !     use bit_rep_data, only: niftot
    !     use SystemData, only: nel
    !
    !     integer(n_int), intent(out) :: ilut(0:niftot)
    !     integer,intent(out) :: nI(nel)
    !     integer :: i
    !
    !     nI = (/ 1,3,6/)
    !     ilut = 0
    !     do i=1,nel
    !         ilut(0) = ibset(ilut(0),nI(i)-1)
    !     enddo
    ! end subroutine generate_test_ilut


end program
