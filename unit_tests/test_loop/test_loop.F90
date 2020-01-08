module test_loop_helpers
    use mpi
    use util_mod, only: get_free_unit
    implicit none
#include "NECICore.h"
    private
    public :: test_loop_factory, FciDumpWriter_t, InputWriter_t, Writer_t

    abstract interface
        subroutine to_unit_writer_t(iunit)
            integer, intent(in) :: iunit
        end subroutine
    end interface

    type, abstract :: Writer_t
        procedure(to_unit_writer_t), pointer, nopass :: write
        character(:), allocatable :: filepath
    end type

    type, extends(Writer_t) :: FciDumpWriter_t
    end type

    type, extends(Writer_t) :: InputWriter_t
    end type

contains

    subroutine test_loop_factory(input_writer, fcidump_writer, additional_writers)
        type(InputWriter_t), intent(in) :: input_writer
        type(FciDumpWriter_t), intent(in) :: fcidump_writer
        class(Writer_t), intent(in), optional :: additional_writers(:)
        integer :: i, myrank, err

        call mpi_comm_rank(MPI_COMM_WORLD, myrank, err)

        if (myrank == 0) then
            call write_file(input_writer)
            call write_file(fcidump_writer)
            if (present(additional_writers)) then
                do i = lbound(additional_writers, 1), ubound(additional_writers, 1)
                    call write_file(additional_writers(i))
                end do
            end if
        end if

        do i = 1, 3
            call NECICore(filename_in=input_writer%filepath, &
                          & int_name=fcidump_writer%filepath, &
                          & call_as_lib=.true.)
        end do

        if (myrank == 0) then
            call delete_file(input_writer%filepath)
            call delete_file(fcidump_writer%filepath)
            if (present(additional_writers)) then
                do i = lbound(additional_writers, 1), ubound(additional_writers, 1)
                    call delete_file(additional_writers(i)%filepath)
                end do
            end if
        end if
    end subroutine test_loop_factory

    subroutine delete_file(path)
        character(*), intent(in) :: path
        integer :: file_id

        ! Even though the file is only deleted, the unit has to be assigned
        file_id = get_free_unit()
        
        open(file_id, file=path, status='old')
        close(file_id, status='delete')
    end subroutine

    subroutine write_file(writer)
        class(Writer_t), intent(in) :: writer
        integer :: file_id

        file_id = get_free_unit()
        open(file_id, file=writer%filepath)
            call writer%write(file_id)
        close(file_id)
    end subroutine

end module test_loop_helpers

module test_loop_testcases
    use test_loop_helpers, only: &
        test_loop_factory, InputWriter_t, FciDumpWriter_t
    implicit none
    private
    public :: test_loop, test_loop_pchb

contains

    !> FCIDUMP file for 1s and 2s orbital of He-atom.
    subroutine write_He_fcidump(unit_id)
        integer, intent(in) :: unit_id
        write(unit_id, '(A)') ' &FCI NORB=2,NELEC=2,MS2=0,'
        write(unit_id, '(A)') '  ORBSYM=1,1,'
        write(unit_id, '(A)') '  ISYM=0'
        write(unit_id, '(A)') ' &END'
        write(unit_id, '(A)') '     1.0232804328        1    1    1    1'
        write(unit_id, '(A)') '   -0.32447343993        2    1    1    1'
        write(unit_id, '(A)') '    0.24653188123        2    1    2    1'
        write(unit_id, '(A)') '    0.89791396824        2    2    1    1'
        write(unit_id, '(A)') '   -0.30190988085        2    2    2    1'
        write(unit_id, '(A)') '    0.84005597186        2    2    2    2'
        write(unit_id, '(A)') '    -1.9424691161        1    1    0    0'
        write(unit_id, '(A)') '    0.32447343993        2    1    0    0'
        write(unit_id, '(A)') '    0.52336359664E-02    2    2    0    0'
        write(unit_id, '(A)') '   -0.91919000000        1    0    0    0'
        write(unit_id, '(A)') '     1.5545000000        2    0    0    0'
        write(unit_id, '(A)') '     0.0000000000        0    0    0    0'
    end subroutine

    subroutine test_loop_pchb()
        call test_loop_factory(InputWriter_t(create_input, 'NECI_input'), &
                               FciDumpWriter_t(write_He_fcidump, 'FCIDUMP'))
        contains

        subroutine create_input(unit_id)
            integer, intent(in) :: unit_id
            write(unit_id, '(A)') 'Title'

            write(unit_id, '(A)') 'System read'
            write(unit_id, '(A)') '    electrons  2'
            write(unit_id, '(A)') '    nonuniformrandexcits pchb'
            write(unit_id, '(A)') '    nobrillouintheorem'
            write(unit_id, '(A)') '    freeformat'
            write(unit_id, '(A)') 'endsys'

            write(unit_id, '(A)') 'calc'

            write(unit_id, '(A)') '    totalwalkers 1000'
            write(unit_id, '(A)') '    (readpops'
            write(unit_id, '(A)') '    (walkcontgrow'
            write(unit_id, '(A)') '    semi-stochastic 10'

            write(unit_id, '(A)') '    methods'
            write(unit_id, '(A)') '        method vertex fcimc'
            write(unit_id, '(A)') '    endmethods'

            write(unit_id, '(A)') '    diagshift .00'
            write(unit_id, '(A)') '    shiftdamp .02'
            write(unit_id, '(A)') '    nmcyc 50000'
            write(unit_id, '(A)') '    stepsshift 10'
            write(unit_id, '(A)') '    proje-changeref 1.20'
            write(unit_id, '(A)') '    truncinitiator'
            write(unit_id, '(A)') '    addtoinitiator  3'
            write(unit_id, '(A)') '    allrealcoeff'
            write(unit_id, '(A)') '    realspawncutoff .30'
            write(unit_id, '(A)') '    jump-shift'
            write(unit_id, '(A)') '    tau 0.01 search'
            write(unit_id, '(A)') '    max-tau .02'
            write(unit_id, '(A)') '    maxwalkerbloom 1'
            write(unit_id, '(A)') '    memoryfacspawn 10.00'
            write(unit_id, '(A)') '    memoryfacpart 5.00'
            write(unit_id, '(A)') '    time 200'
            write(unit_id, '(A)') '    startsinglepart 10'
            write(unit_id, '(A)') '    pops-core 10000'
            write(unit_id, '(A)') '    rdmsamplingiters 30'
            write(unit_id, '(A)') 'endcalc'

            write(unit_id, '(A)') 'logging'
            write(unit_id, '(A)') '    highlypopwrite 50'
            write(unit_id, '(A)') '    print-spin-resolved-RDMs'
            write(unit_id, '(A)') '    printonerdm'
            write(unit_id, '(A)') '    calcrdmonfly 3 10 10'
            write(unit_id, '(A)') 'endlog'
            write(unit_id, '(A)') 'end'
        end subroutine create_input
    end subroutine test_loop_pchb

    subroutine test_loop()
        call test_loop_factory(InputWriter_t(create_input, 'NECI_input'), &
                               FciDumpWriter_t(write_He_fcidump, 'FCIDUMP'))
        contains

        subroutine create_input(unit_id)
            integer, intent(in) :: unit_id
            write(unit_id, '(A)') 'Title'

            write(unit_id, '(A)') 'System read'
            write(unit_id, '(A)') '    electrons  2'
            write(unit_id, '(A)') '    nonuniformrandexcits 4ind-weighted-2'
            write(unit_id, '(A)') '    nobrillouintheorem'
            write(unit_id, '(A)') '    freeformat'
            write(unit_id, '(A)') 'endsys'

            write(unit_id, '(A)') 'calc'

            write(unit_id, '(A)') '    totalwalkers 1000'
            write(unit_id, '(A)') '    (readpops'
            write(unit_id, '(A)') '    (walkcontgrow'
            write(unit_id, '(A)') '    semi-stochastic 10'

            write(unit_id, '(A)') '    methods'
            write(unit_id, '(A)') '        method vertex fcimc'
            write(unit_id, '(A)') '    endmethods'

            write(unit_id, '(A)') '    diagshift .00'
            write(unit_id, '(A)') '    shiftdamp .02'
            write(unit_id, '(A)') '    nmcyc 50000'
            write(unit_id, '(A)') '    stepsshift 10'
            write(unit_id, '(A)') '    proje-changeref 1.20'
            write(unit_id, '(A)') '    truncinitiator'
            write(unit_id, '(A)') '    addtoinitiator  3'
            write(unit_id, '(A)') '    allrealcoeff'
            write(unit_id, '(A)') '    realspawncutoff .30'
            write(unit_id, '(A)') '    jump-shift'
            write(unit_id, '(A)') '    tau 0.01 search'
            write(unit_id, '(A)') '    max-tau .02'
            write(unit_id, '(A)') '    maxwalkerbloom 1'
            write(unit_id, '(A)') '    memoryfacspawn 10.00'
            write(unit_id, '(A)') '    memoryfacpart 5.00'
            write(unit_id, '(A)') '    time 200'
            write(unit_id, '(A)') '    startsinglepart 10'
            write(unit_id, '(A)') '    pops-core 10000'
            write(unit_id, '(A)') '    rdmsamplingiters 30'
            write(unit_id, '(A)') 'endcalc'

            write(unit_id, '(A)') 'logging'
            write(unit_id, '(A)') '    highlypopwrite 50'
            write(unit_id, '(A)') '    print-spin-resolved-RDMs'
            write(unit_id, '(A)') '    printonerdm'
            write(unit_id, '(A)') '    calcrdmonfly 3 10 10'
            write(unit_id, '(A)') 'endlog'
            write(unit_id, '(A)') 'end'
        end subroutine create_input
    end subroutine test_loop
end module test_loop_testcases

program test_loop_program

    use mpi
    use fruit
    use test_loop_testcases, only: test_loop, test_loop_pchb

    implicit none
    integer :: failed_count, err

    call mpi_init(err)

    call init_fruit()

    call test_loop_driver()

    call fruit_summary()
    call fruit_finalize()
    call get_failed_count(failed_count)

    if (failed_count /= 0) call stop_all('test_loop_program', 'failed_tests')

    call mpi_finalize(err)

contains

    subroutine test_loop_driver()
        call run_test_case(test_loop, "test_loop")
        call run_test_case(test_loop_pchb, "test_loop_pchb")
    end subroutine test_loop_driver
end program test_loop_program
