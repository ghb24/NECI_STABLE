#include "macros.h"

module block_spawns

    use SystemData, only: nel
    use constants, only: n_int, stdout
    use util_mod, only: stop_all, get_free_unit
    use input_parser_mod, only: TokenIterator_t, FileReader_t, parse_definedet,&
                                ManagingFileReader_t
    use DetBitOps, only: EncodeBitDet
    use bit_rep_data, only: nifd
    use bit_reps, only: decode_bit_det
    use FciMCData, only: ll_node
    use hash, only: init_hash_table, fill_in_hash_table, clear_hash_table, &
                    hash_table_lookup

    better_implicit_none

    private

    public :: init_block_list, tBlockSpawns, blocklist_name, block_list, &
              is_excitation_aborted, blocklist_hash_table

    logical :: tBlockSpawns
    character(len=255) :: blocklist_name
    integer(kind=n_int), allocatable :: block_list(:, :)
    type(ll_node), pointer :: blocklist_hash_table(:)

contains

    subroutine init_block_list(blocklistname, blocklist, blocklist_ht)
        integer :: iunit, numlines, i, stat, nI_blocked(nel)
        integer(kind=n_int), allocatable, intent(inout) :: blocklist(:,:)
        type(ll_node), pointer, intent(inout) :: blocklist_ht(:)
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
        allocate(blocklist(0:nifd, numlines), source=0_n_int)

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

        call init_blocklist_hashtable(blocklist, blocklist_ht)

    end subroutine init_block_list

    subroutine init_blocklist_hashtable(blocklist, blocklist_ht)
        integer(kind=n_int), intent(in) :: blocklist(0:,:)
        type(ll_node), pointer, intent(inout) :: blocklist_ht(:)

            allocate(blocklist_ht(size(blocklist, 2)))
            call init_hash_table(blocklist_ht)
            call fill_in_hash_table(blocklist_ht, size(blocklist, 2), &
                                    blocklist, size(blocklist, 2), .false.)

    end subroutine init_blocklist_hashtable

    function is_excitation_aborted(excited_nJ, excited_ilut) result(tAbort)
        integer, intent(in) :: excited_nJ(nel)
        integer(kind=n_int), intent(in) :: excited_ilut(0:nifd)
        integer(kind=n_int), allocatable :: blocklist(:,:)
        logical :: tAbort
        integer :: i, blocked_nJ(nel), idx, DetHash

        tAbort = .false.

        if (allocated(block_list)) then
            call hash_table_lookup(excited_nJ, excited_ilut, nifd, &
                                   blocklist_hash_table, block_list, idx, &
                                   DetHash, tAbort)
        end if

    end function is_excitation_aborted

end module block_spawns
