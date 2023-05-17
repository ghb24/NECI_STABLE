#include "macros.h"

module test_loop_testcases
    use test_loop_helpers, only: test_loop_factory
    use unit_test_helper_excitgen, only: FciDumpWriter_t, InputWriter_t
    implicit none
    private
    public :: test_loop_4ind_wghtd_2, test_loop_pchb, test_loop_guga

    abstract interface
        subroutine ptr_to_unit_writer_t(iunit)
            integer, intent(in) :: iunit
        end subroutine
    end interface


    type, extends(InputWriter_t) :: RunTimePointerInputWriter_t
        procedure(ptr_to_unit_writer_t), nopass, pointer :: ptr_write_to_unit => null()
    contains
        procedure :: write_to_unit
    end type

    type, extends(FciDumpWriter_t) :: HeFciDumpWriter_t
    contains
        procedure :: write_to_unit => HeFciDumpWriter_t_write_to_unit
    end type

contains

    ! In the long run, this procedure can and should be generalized
    ! if more tests are performed.
    ! Have a look at rasscf::fciqcm_make_inp in the OpenMolcas codebase
    ! for inspiration.
    ! For the time being I stay with YAGNI because I don't know for sure
    ! if we need more tests.
    subroutine create_input(unit_id, exc_generator, flag_guga)
        integer, intent(in) :: unit_id
        character(*), intent(in) :: exc_generator
        logical, intent(in), optional :: flag_guga

        logical :: flag_guga_

        def_default(flag_guga_,flag_guga,.false.)

        write(unit_id, '(A)') 'Title'

        write(unit_id, '(A)') 'System read'
        write(unit_id, '(A)') '    electrons  2'
        write(unit_id, '(A)') '    nonuniformrandexcits '//exc_generator
        write(unit_id, '(A)') '    nobrillouintheorem'
        if (flag_guga_) then
            write(unit_id, '(A)') '    guga 0'
        end if
        write(unit_id, '(A)') '    freeformat'
        write(unit_id, '(A)') 'endsys'

        write(unit_id, '(A)') 'calc'

        write(unit_id, '(A)') '    totalwalkers 100'
        write(unit_id, '(A)') '    # readpops'
        write(unit_id, '(A)') '    # walkcontgrow'
        write(unit_id, '(A)') '    semi-stochastic 10'

        write(unit_id, '(A)') '    methods'
        write(unit_id, '(A)') '        method vertex fcimc'
        write(unit_id, '(A)') '    endmethods'

        write(unit_id, '(A)') '    diagshift .00'
        write(unit_id, '(A)') '    shiftdamp .02'
        write(unit_id, '(A)') '    nmcyc 10000'
        write(unit_id, '(A)') '    stepsshift 10'
        write(unit_id, '(A)') '    proje-changeref 1.20'
        write(unit_id, '(A)') '    truncinitiator'
        write(unit_id, '(A)') '    addtoinitiator  3'
        write(unit_id, '(A)') '    allrealcoeff'
        write(unit_id, '(A)') '    realspawncutoff .30'
        write(unit_id, '(A)') '    jump-shift'
        write(unit_id, '(A)') '    tau-values start tau-factor 0.05'
        write(unit_id, '(A)') '    tau-search algorithm conventional maxwalkerbloom 1.0'
        write(unit_id, '(A)') '    memoryfacspawn 10.00'
        write(unit_id, '(A)') '    memoryfacpart 5.00'
        write(unit_id, '(A)') '    time 200'
        write(unit_id, '(A)') '    startsinglepart 10'
        write(unit_id, '(A)') '    pops-core 10'
        if (.not. flag_guga_) then
            write(unit_id, '(A)') '    rdmsamplingiters 30'
        end if
        write(unit_id, '(A)') 'endcalc'

        write(unit_id, '(A)') 'logging'
        write(unit_id, '(A)') '    highlypopwrite 50'
        if (.not. flag_guga_) then
            write(unit_id, '(A)') '    biased-RDMs'
            write(unit_id, '(A)') '    print-spin-resolved-RDMs'
            write(unit_id, '(A)') '    printonerdm'
            write(unit_id, '(A)') '    calcrdmonfly 3 10 10'
        end if
        write(unit_id, '(A)') 'endlog'
        write(unit_id, '(A)') 'end'
    end subroutine create_input


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
        call test_loop_factory(&
                RunTimePointerInputWriter_t(ptr_write_to_unit=create_input_pchb, filepath='NECI_input'), &
                HeFciDumpWriter_t(filepath='FCIDUMP'))

        contains

        subroutine create_input_pchb(unit_id)
            integer, intent(in) :: unit_id
            call create_input(unit_id, exc_generator='pchb')
        end subroutine
    end subroutine test_loop_pchb

    subroutine test_loop_4ind_wghtd_2()
        call test_loop_factory(&
                RunTimePointerInputWriter_t(filepath='NECI_input', ptr_write_to_unit=create_input_4ind_wghtd_2), &
                HeFciDumpWriter_t(filepath='FCIDUMP'))

        contains

        subroutine create_input_4ind_wghtd_2(unit_id)
            integer, intent(in) :: unit_id
            call create_input(unit_id, exc_generator='4ind-weighted-2')
        end subroutine
    end subroutine test_loop_4ind_wghtd_2

    subroutine test_loop_guga()
        call test_loop_factory(&
            RunTimePointerInputWriter_t(filepath='NECI_input', ptr_write_to_unit=create_input_guga), &
            HeFciDumpWriter_t(filepath='FCIDUMP'))

        contains

            subroutine create_input_guga(unit_id)
                integer, intent(in) :: unit_id
                call create_input(unit_id, exc_generator = 'mol_guga_weighted', &
                                  flag_guga = .true.)
            end subroutine create_input_guga
    end subroutine test_loop_guga

    subroutine write_to_unit(this, iunit)
        class(RunTimePointerInputWriter_t), intent(in) :: this
        integer, intent(in) :: iunit
        call this%ptr_write_to_unit(iunit)
    end subroutine

    subroutine HeFciDumpWriter_t_write_to_unit(this, iunit)
        class(HeFciDumpWriter_t), intent(in) :: this
        integer, intent(in) :: iunit
        unused_var(this)
        call write_He_fcidump(iunit)
    end subroutine
end module test_loop_testcases

