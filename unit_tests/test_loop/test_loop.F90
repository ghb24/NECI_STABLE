program test_loop_program

    use mpi
    use fruit

    use util_mod, only: get_free_unit

    implicit none
#include "NECICore.h"

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
    end subroutine test_loop_driver

    subroutine test_loop()
        character(*), parameter :: fcidump = 'FCIDUMP', input = 'NECI_input'
        integer :: fcidump_id, input_id
        integer :: i, myrank, err

        call mpi_comm_rank(MPI_COMM_WORLD, myrank, err)

        if (myrank == 0) then
            input_id = get_free_unit()
            open(input_id, file=input, status='new')
            call create_input(input_id)
            close(input_id)

            fcidump_id = get_free_unit()
            open(fcidump_id, file=fcidump, status='new')
            call create_fcidump(fcidump_id)
            close(fcidump_id)
        end if

        do i = 1, 3
            call NECICore(filename_in=input, int_name=fcidump, &
                          & call_as_lib=.true.)
        end do

        if (myrank == 0) then
            open(input_id, file=input, status='old')
            close(input_id, status='delete')
            open(fcidump_id, file=fcidump, status='old')
            close(fcidump_id, status='delete')
        end if

    end subroutine test_loop

    subroutine create_fcidump(unit_id)
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
    end subroutine

end program test_loop_program
