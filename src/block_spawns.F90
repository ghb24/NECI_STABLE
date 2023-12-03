#include "macros.h"

module block_spawns

    use SystemData, only: nel
    use constants, only: n_int, stdout
    use util_mod, only: stop_all, get_free_unit
    use input_parser_mod, only: TokenIterator_t, FileReader_t, parse_definedet,&
                                ManagingFileReader_t
    use DetBitOps, only: EncodeBitDet
    use bit_rep_data, only: niftot
    use bit_reps, only: decode_bit_det

    better_implicit_none

    private

    public :: init_block_list, tBlockSpawns, blocklist_name, block_list, &
              is_excitation_aborted

    logical :: tBlockSpawns
    character(len=255) :: blocklist_name
    integer(kind=n_int), allocatable :: block_list(:, :)

contains

    subroutine init_block_list(blocklistname, blocklist)
        integer :: iunit, numlines, i, stat, nI_blocked(nel)
        integer(kind=n_int), allocatable :: blocklist(:,:)
        character(len=255), intent(in) :: blocklistname
        character(len=500) :: line
        type(TokenIterator_t) :: tokens
        character(*), parameter :: t_r = 'init_block_list'
        class(FileReader_t), allocatable :: filereader

! has to be replaced with growing buffer
        ! Count the number of lines
        iunit = get_free_unit()
        open(newunit=iunit, file=blocklistname)
        numlines = 0
        do
            read(iunit, *, iostat=stat) line
            if (stat /= 0) exit
            numlines = numlines + 1
        end do
        close(iunit)
        if (numlines == 0) then
            call stop_all(t_r, trim(blocklistname)//" file is empty.")
        end if
        allocate(blocklist(0:niftot, numlines), source=0_n_int)

        filereader = ManagingFileReader_t(blocklistname)
        do i = 1, numlines
            if (filereader%nextline(tokens, skip_empty=.false.)) then
                call parse_definedet(tokens, nI_blocked)
                call EncodeBitDet(nI_blocked, blocklist(:, i))
            else
                call stop_all(t_r, "Unexpected EOF reached.")
            end if
#ifdef DEBUG_
            write(stdout, *) "Determinant/CSF", i, "to be blocked:", nI_blocked
#endif
        end do
        call filereader%close()

        write(stdout, *) 'Read in', numlines, 'determinants/CSFs to be blocked.'

    end subroutine init_block_list

    function is_excitation_aborted(excited_nJ, blocklist) result(tAbort)
        integer, intent(in) :: excited_nJ(nel)
        integer(kind=n_int), intent(in), optional, allocatable :: blocklist(:,:)
        logical :: tAbort
        integer :: i, blocked_nJ(nel)

        tAbort = .false.

        if (present(blocklist)) then
            do i = 1, size(blocklist, 2)
                call decode_bit_det(blocked_nJ, blocklist(:, i))
                if (all(blocked_nJ == excited_nJ)) then
                    tAbort = .true.
                    exit
                end if
            end do
        end if

    end function is_excitation_aborted

end module block_spawns